#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/Format.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/misc/ExactStencilWeights.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/GaussSeidelSmoother.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "../../data/param/gmg_P2.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t minLevel         = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel         = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t max_outer_iter   = parameters.getParameter< uint_t >( "max_outer_iter" );
   //const uint_t max_cg_iter      = parameters.getParameter< uint_t >( "max_cg_iter" );
   const real_t mg_tolerance     = parameters.getParameter< real_t >( "mg_tolerance" );
   //const real_t coarse_tolerance = parameters.getParameter< real_t >( "coarse_tolerance" );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( parameters.getParameter< std::string >( "mesh" ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   hhg::P2Function< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P2Function< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P2Function< real_t > Lu( "Lu", storage, minLevel, maxLevel );
   hhg::P2Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   hhg::P2Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hhg::P2Function< real_t > err( "err", storage, minLevel, maxLevel );
   hhg::P2Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hhg::Point3D& ) > rhs   = []( const hhg::Point3D& ) { return 0; };
   std::function< real_t( const hhg::Point3D& ) > zero  = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones  = []( const hhg::Point3D& ) { return 1.0; };
   walberla::math::seedRandomGenerator( 0 );
   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) { return walberla::math::realRandom( 0.0, 20.0 ); };

   WALBERLA_LOG_INFO_ON_ROOT( "Interpolating u" );
   u.interpolate( rand, maxLevel, hhg::Inner );
   u.interpolate( exact, maxLevel, hhg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );

   //  WALBERLA_LOG_INFO_ON_ROOT("Interpolating and integrating rhs");
   //  npoints_helper.interpolate(rhs, maxLevel);
   //  M.apply(npoints_helper, f, maxLevel, hhg::All);

   WALBERLA_LOG_INFO_ON_ROOT( "Setting up stiffness operator" );
   auto                           start = walberla::timing::getWcTime();
   hhg::P2ConstantLaplaceOperator L( storage, minLevel, maxLevel );
   auto                           end       = walberla::timing::getWcTime();
   real_t                         setupTime = end - start;

   npoints_helper.interpolate( ones, maxLevel );
   real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

   auto smoother = std::make_shared< hhg::GaussSeidelSmoother<hhg::P2ConstantLaplaceOperator>  >();
   auto coarseGridSolver = std::make_shared< hhg::CGSolver< hhg::P2ConstantLaplaceOperator > >( storage, minLevel, minLevel );
   auto restrictionOperator = std::make_shared< hhg::P2toP2QuadraticRestriction>();
   auto prolongationOperator = std::make_shared< hhg::P2toP2QuadraticProlongation >();

   auto gmgSolver = hhg::GeometricMultigridSolver< hhg::P2ConstantLaplaceOperator >(
      storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

   if( parameters.getParameter< bool >( "useExactWeights" ) )
   {
      WALBERLA_LOG_INFO( "WARNING: works only on tri_1el mesh" );
      auto weights = hhg::stencilWeights::tri_1el();

      for( uint_t i = minLevel; i <= maxLevel; ++i )
      {
         real_t* vToV =
             storage->getFace( PrimitiveID( 6 ) )->getData( L.getVertexToVertexOpr().getFaceStencilID() )->getPointer( i );
         for( uint_t j = 0; j < 7; ++j )
         {
            vToV[j] = weights.vertexToVertexStencil[j];
         }

         real_t* eToV =
             storage->getFace( PrimitiveID( 6 ) )->getData( L.getEdgeToVertexOpr().getFaceStencilID() )->getPointer( i );
         for( uint_t j = 0; j < 12; ++j )
         {
            eToV[j] = weights.edgeToVertexStencil[j];
         }

         real_t* vToE =
             storage->getFace( PrimitiveID( 6 ) )->getData( L.getVertexToEdgeOpr().getFaceStencilID() )->getPointer( i );
         for( uint_t j = 0; j < 12; ++j )
         {
            vToE[j] = weights.vertexToEdgeStencil[j];
         }

         real_t* eToE = storage->getFace( PrimitiveID( 6 ) )->getData( L.getEdgeToEdgeOpr().getFaceStencilID() )->getPointer( i );
         for( uint_t j = 0; j < 15; ++j )
         {
            eToE[j] = weights.edgeToEdgeStencil[j];
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Starting V cycles" );
   WALBERLA_LOG_INFO_ON_ROOT(
       hhg::format( "%6s|%10s|%10s|%10s|%10s|%10s", "iter", "abs_res", "rel_res", "conv", "L2-error", "Time" ) );

   real_t rel_res = 1.0;

   L.apply( u, Lu, maxLevel, hhg::Inner );
   r.assign( {1.0, -1.0}, {f, Lu}, maxLevel, hhg::Inner );

   real_t begin_res   = std::sqrt( r.dotGlobal( r, maxLevel, hhg::Inner ) );
   real_t abs_res_old = begin_res;

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

   //WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  -", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err));
   WALBERLA_LOG_INFO_ON_ROOT(
       hhg::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res / abs_res_old, discr_l2_err, 0 ) )

   real_t       solveTime              = real_c( 0.0 );
   real_t       averageConvergenceRate = real_c( 0.0 );
   const uint_t convergenceStartIter   = 3;

   uint_t i = 0;
   for( ; i < max_outer_iter; ++i )
   {
      start = walberla::timing::getWcTime();

      // Apply P2 geometric multigrid solver
      gmgSolver.solve( L, u, f, maxLevel );

      end = walberla::timing::getWcTime();

      L.apply( u, Lu, maxLevel, hhg::Inner );
      r.assign( {1.0, -1.0}, {f, Lu}, maxLevel, hhg::Inner );
      real_t abs_res = std::sqrt( r.dotGlobal( r, maxLevel, hhg::Inner ) );
      rel_res        = abs_res / begin_res;
      err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
      discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

      //WALBERLA_LOG_INFO_ON_ROOT(fmt::format("{:3d}   {:e}  {:e}  {:e}  {:e}  {:e}", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err, end-start));
      WALBERLA_LOG_INFO_ON_ROOT( hhg::format(
          "%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", i + 1, abs_res, rel_res, abs_res / abs_res_old, discr_l2_err, end - start ) )
      solveTime += end - start;

      if( i >= convergenceStartIter )
      {
         averageConvergenceRate += abs_res / abs_res_old;
      }

      abs_res_old = abs_res;

      if( rel_res < mg_tolerance )
      {
         break;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "MinLevel: " << minLevel );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "MaxLevel: " << maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "DoFs on MaxLevel: " << (uint_t) npoints );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "Setup time: " << std::defaultfloat << setupTime );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "Solve time: " << std::defaultfloat << solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "Time to solution: " << std::defaultfloat << setupTime + solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "Avg. convergence rate: " << std::scientific
                                              << averageConvergenceRate / real_c( i - convergenceStartIter ) );
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 25 ) << "L^2 error: " << std::scientific << discr_l2_err );

   WALBERLA_CHECK_LESS( discr_l2_err, 8.4e-08 )

   if( parameters.getParameter< bool >( "vtkOutput" ) )
   {
      VTKOutput vtkOutput("../output", "gmg_P2_h_refinement", storage);
      vtkOutput.add( u );
      vtkOutput.add( u_exact );
      vtkOutput.add( f );
      vtkOutput.add( r );
      vtkOutput.add( err );
      vtkOutput.add( npoints_helper );
      vtkOutput.write( maxLevel );
   }

   if( parameters.getParameter< bool >( "printTiming" ) )
   {
      walberla::WcTimingTree tt = timingTree->getReduced();
      WALBERLA_LOG_INFO_ON_ROOT( tt );
   }

   return 0;
}
