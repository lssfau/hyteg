#include "core/config/Create.h"
#include "core/DataTypes.h"

#include "tinyhhg_core/composites/P1BlendingStokesOperator.hpp"
#include "tinyhhg_core/composites/P1CoefficientStokesOperator.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/Format.hpp"
#include "tinyhhg_core/geometry/CircularMap.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1BlendingOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSmoother.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"

using walberla::real_t;
using namespace hhg;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg;

   if( argc == 1 )
   {
      walberla::shared_ptr< walberla::config::Config > cfg_( new walberla::config::Config );
      cfg_->readParameterFile( "../../data/param/stokes_blending.prm" );
      cfg = cfg_;
   } else
   {
      cfg = walberla::config::create( argc, argv );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "config = " << *cfg );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const std::string meshFileName  = parameters.getParameter< std::string >( "meshFileName" );
   const uint_t      minLevel      = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel      = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t      coarseMaxiter = parameters.getParameter< uint_t >( "coarseMaxiter" );
   const uint_t      maxOuterIter  = parameters.getParameter< uint_t >( "maxOuterIter" );

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{
       {parameters.getParameter< real_t >( "circleCenterX" ), parameters.getParameter< real_t >( "circleCenterY" ), 0}};
   real_t circleRadius = parameters.getParameter< real_t >( "circleRadius" );

   for( auto it : setupStorage.getFaces() )
   {
      Face& face = *(it.second);

      std::vector< PrimitiveID > neighborEdgesOnBoundary = face.neighborEdges();
      std::remove_if( neighborEdgesOnBoundary.begin(), neighborEdgesOnBoundary.end(), [&setupStorage]( const PrimitiveID& id ) {
         return !setupStorage.onBoundary( id );
      } );

      if( neighborEdgesOnBoundary.size() > 0 )
      {
         Edge& edge = *setupStorage.getEdge( neighborEdgesOnBoundary[0] );

         if( ( edge.getCoordinates()[0] - circleCenter ).norm() < 0.4 )
         {
            setupStorage.setGeometryMap( edge.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
            setupStorage.setGeometryMap( face.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
         }
      }
   }

   hhg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::P1Function< real_t >       one( "one", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > tmp2( "tmp", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > err( "err", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > Lu( "Lu", storage, minLevel, maxLevel );

   auto smoother = std::make_shared< hhg::UzawaSmoother< hhg::P1BlendingStokesOperator > >(
       storage, minLevel, maxLevel, 0.3 );
   auto coarseGridSolver = std::make_shared< hhg::MinResSolver< hhg::P1BlendingStokesOperator > >( storage, minLevel, minLevel, coarseMaxiter );
   auto restrictionOperator = std::make_shared< hhg::P1P1StokesToP1P1StokesRestriction>();
   auto prolongationOperator = std::make_shared< hhg::P1P1StokesToP1P1StokesProlongation >();

   auto solver = hhg::GeometricMultigridSolver< hhg::P1BlendingStokesOperator >(
      storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2, 2 );

   auto          start = walberla::timing::getWcTime();
   hhg::P1BlendingStokesOperator L( storage, minLevel, maxLevel );
   auto          end       = walberla::timing::getWcTime();
   real_t        setupTime = end - start;

   P1BlendingMassOperator M( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > zeros   = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones    = []( const hhg::Point3D& ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > exact_u = []( const hhg::Point3D& x ) {
      return sin( x[0] ) * cos( x[1] ) / ( x[0] + 1 );
   };
   std::function< real_t( const hhg::Point3D& ) > exact_v = []( const hhg::Point3D& x ) {
      return -( ( x[0] + 1 ) * cos( x[0] ) - sin( x[0] ) ) * sin( x[1] ) / pow( x[0] + 1, 2 );
   };
   std::function< real_t( const hhg::Point3D& ) > exact_p = []( const hhg::Point3D& x ) {
      return pow( x[0], 2 ) * pow( x[1], 3 );
   };
   std::function< real_t( const hhg::Point3D& ) > rhs_u = []( const hhg::Point3D& x ) {
      return ( 2 * x[0] * pow( x[1], 3 ) * pow( x[0] + 1, 3 ) + 2.0 * pow( x[0] + 1, 2 ) * sin( x[0] ) * cos( x[1] ) +
               ( 4.0 * x[0] + 4.0 ) * cos( x[0] ) * cos( x[1] ) -
               2.0 * ( ( x[0] + 1 ) * cos( x[0] ) + sin( x[0] ) ) * cos( x[1] ) ) /
             pow( x[0] + 1, 3 );
   };
   std::function< real_t( const hhg::Point3D& ) > rhs_v = []( const hhg::Point3D& x ) {
      return ( 3 * pow( x[0], 2 ) * pow( x[1], 2 ) * pow( x[0] + 1, 4 ) +
               2.0 * pow( x[0] + 1, 2 ) * ( -( x[0] + 1 ) * cos( x[0] ) + 2 * sin( x[0] ) ) * sin( x[1] ) +
               6.0 * ( ( x[0] + 1 ) * cos( x[0] ) - sin( x[0] ) ) * sin( x[1] ) ) /
             pow( x[0] + 1, 4 );
   };

   u.u.interpolate( exact_u, maxLevel, hhg::DirichletBoundary );
   u.v.interpolate( exact_v, maxLevel, hhg::DirichletBoundary );

   hhg::communication::syncFunctionBetweenPrimitives( u.u, maxLevel );
   hhg::communication::syncFunctionBetweenPrimitives( u.v, maxLevel );

   u_exact.u.interpolate( exact_u, maxLevel, hhg::All );
   u_exact.v.interpolate( exact_v, maxLevel, hhg::All );
   u_exact.p.interpolate( exact_p, maxLevel, hhg::All );

   hhg::communication::syncFunctionBetweenPrimitives( u_exact.u, maxLevel );
   hhg::communication::syncFunctionBetweenPrimitives( u_exact.v, maxLevel );
   hhg::communication::syncFunctionBetweenPrimitives( u_exact.p, maxLevel );

   hhg::vertexdof::projectMean( u_exact.p, maxLevel );

   // Integrate RHS for u
   tmp.u.interpolate( rhs_u, maxLevel, hhg::All );
   //hhg::communication::syncFunctionBetweenPrimitives( tmp.u, maxLevel );
   M.apply( tmp.u, f.u, maxLevel, hhg::All );

   // Integrate RHS for v
   tmp.v.interpolate( rhs_v, maxLevel, hhg::All );
   //hhg::communication::syncFunctionBetweenPrimitives( tmp.v, maxLevel );
   M.apply( tmp.v, f.v, maxLevel, hhg::All );

   // Apply compatibility projection
   one.interpolate( zeros, maxLevel, hhg::All );
   one.interpolate( ones, maxLevel, hhg::DirichletBoundary );
   communication::syncFunctionBetweenPrimitives( one, maxLevel);
   L.div_x.apply( u_exact.u, tmp2.p, maxLevel, hhg::DirichletBoundary, Replace );
   L.div_y.apply( u_exact.v, tmp2.p, maxLevel, hhg::DirichletBoundary, Add );
   real_t corr = one.dotGlobal( tmp2.p, maxLevel, hhg::DirichletBoundary );
   M.apply( one, tmp2.p, maxLevel, hhg::DirichletBoundary, Replace );
   real_t volume = one.dotGlobal( tmp2.p, maxLevel, hhg::DirichletBoundary );
   tmp2.p.assign( {corr / volume}, {one}, maxLevel, hhg::DirichletBoundary );
   M.apply( tmp2.p, f.p, maxLevel, hhg::DirichletBoundary );

   one.interpolate( ones, maxLevel, hhg::All );
   communication::syncFunctionBetweenPrimitives( one, maxLevel);
   real_t npoints = one.dotGlobal( one, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Starting Uzawa cycles" );
   WALBERLA_LOG_INFO_ON_ROOT( hhg::format( "%6s|%10s|%10s|%10s|%10s", "iter", "abs_res", "rel_res", "conv", "Time" ) );

   L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, r}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   real_t begin_res   = std::sqrt( r.dotGlobal( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) );
   real_t abs_res_old = begin_res;
   real_t rel_res     = 1.0;

   WALBERLA_LOG_INFO_ON_ROOT(
       hhg::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res / abs_res_old, 0 ) );

   real_t       solveTime              = real_c( 0.0 );
   real_t       averageConvergenceRate = real_c( 0.0 );
   const uint_t convergenceStartIter   = 3;

   uint_t outer;
   for( outer = 0; outer < maxOuterIter; ++outer )
   {
      start = walberla::timing::getWcTime();

      solver.solve(L, u, f, maxLevel);

      end = walberla::timing::getWcTime();
      hhg::vertexdof::projectMean( u.p, maxLevel );

      L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );

      r.assign( {1.0, -1.0}, {f, r}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      real_t abs_res = std::sqrt( r.dotGlobal( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) );
      rel_res        = abs_res / begin_res;
      WALBERLA_LOG_INFO_ON_ROOT(
          hhg::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e", outer + 1, abs_res, rel_res, abs_res / abs_res_old, end - start ) );
      solveTime += end - start;

      if( abs_res / abs_res_old > 0.95 )
      {
         break;
      }

      if( outer >= convergenceStartIter )
      {
         averageConvergenceRate += abs_res / abs_res_old;
      }

      abs_res_old = abs_res;

      //      if (rel_res < mg_tolerance)
      //      {
      //         break;
      //      }
   }

   //   WALBERLA_CHECK_LESS( outer, maxOuterIter );

   WALBERLA_LOG_INFO_ON_ROOT( "Setup time: " << setupTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Solve time: " << solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Time to solution: " << setupTime + solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Avg. convergence rate: " << std::scientific
                                                        << averageConvergenceRate / real_c( outer + 1 - convergenceStartIter ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Dofs: " << 3 * npoints );

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
   real_t discr_u_l2_err = std::sqrt( ( err.u.dotGlobal( err.u, maxLevel ) + err.v.dotGlobal( err.v, maxLevel ) ) / ( 2 * npoints ) );
   real_t discr_p_l2_err = std::sqrt( ( err.p.dotGlobal( err.p, maxLevel ) ) / ( npoints ) );

   WALBERLA_LOG_INFO_ON_ROOT( "velocity_err = " << std::scientific << discr_u_l2_err );
   WALBERLA_LOG_INFO_ON_ROOT( "pressure_err = " << std::scientific << discr_p_l2_err );

   WALBERLA_CHECK_LESS( discr_u_l2_err, 2e-4 );
   WALBERLA_CHECK_LESS( discr_p_l2_err, 2.05e-2 );

   // u_u*iHat + u_v*jHat
   //   hhg::VTKOutput vtkOutput( "../output", "stokes_stab_varcoeff" );
   //   vtkOutput.add( u.u );
   //   vtkOutput.add( u.v );
   //   vtkOutput.add( u.p );
   //   vtkOutput.add( u_exact.u );
   //   vtkOutput.add( u_exact.v );
   //   vtkOutput.add( u_exact.p );
   //   vtkOutput.add( err.u );
   //   vtkOutput.add( err.v );
   //   vtkOutput.add( err.p );
   //   vtkOutput.add( r.u );
   //   vtkOutput.add( r.v );
   //   vtkOutput.add( r.p );
   //   vtkOutput.add( f.u );
   //   vtkOutput.add( f.v );
   //   vtkOutput.add( f.p );
   //   vtkOutput.write( maxLevel, 0 );
   return EXIT_SUCCESS;
}
