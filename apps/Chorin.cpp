#include <boost/core/null_deleter.hpp>


#include "tinyhhg_core/dgfunctionspace/DGUpwindOperator.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultiGrid.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   timingTree->start( "Global" );
   std::string meshFileName = "../data/meshes/flow_around_cylinder.msh";

   real_t viscosity = 1e-4;

   bool   neumann  = true;
   uint_t minLevel = 2;
   uint_t maxLevel = 3;

   real_t time              = 0.0;
   real_t inflowBuildupTime = 0.0;
   real_t endTime           = 5.0;
   uint_t iter              = 0;
   uint_t max_cg_iter       = 50;
   uint_t outerIterations   = 2;

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage, timingTree );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
   loadbalancing::distributed::parmetis( *storage );
#endif

   //  const real_t minimalEdgeLength = hhg::MeshQuality::getMinimalEdgeLength(storage, maxLevel);
   //  real_t dt = 0.025 * minimalEdgeLength;
   real_t dt      = 5e-6;
   real_t dt_plot = 0.005;

   uint_t plotModulo = uint_c( std::ceil( dt_plot / dt ) );

   WALBERLA_LOG_INFO_ON_ROOT( "dt = " << dt );

   std::function< real_t( const hhg::Point3D& ) > bc_x = [&time, &inflowBuildupTime]( const hhg::Point3D& x ) {
      const real_t U_m = 5.0;

      if( x[0] < 1e-8 )
      {
         real_t velocity = 4.0 * U_m * x[1] * ( 0.41 - x[1] ) / ( 0.41 * 0.41 );
         real_t damping;

         if( time < inflowBuildupTime )
         {
            damping = 0.5 * ( 1.0 + std::cos( walberla::math::PI * ( time / inflowBuildupTime - 1.0 ) ) );
         } else
         {
            damping = 1.0;
         }

         return damping * velocity;
      } else
      {
         return 0.0;
      }
   };

   std::function< real_t( const hhg::Point3D& ) > bc_y = []( const hhg::Point3D& ) { return 0.0; };

   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > one  = []( const hhg::Point3D& ) { return 1.0; };

   hhg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > v( "v", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > p( "p", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > p_rhs( "p_rhs", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > p_res( "p_res", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > tmp2( "tmp2", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > res( "res", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > ones( "ones", storage, minLevel, maxLevel );

   auto u_dg = std::make_shared< hhg::DGFunction< real_t > >( "u_dg", storage, minLevel, maxLevel );
   auto v_dg = std::make_shared< hhg::DGFunction< real_t > >( "v_dg", storage, minLevel, maxLevel );

   auto u_dg_old = std::make_shared< hhg::DGFunction< real_t > >( "u_dg", storage, minLevel, maxLevel );
   auto v_dg_old = std::make_shared< hhg::DGFunction< real_t > >( "v_dg", storage, minLevel, maxLevel );

   hhg::P1LaplaceOperator A( storage, minLevel, maxLevel );
   hhg::P1LaplaceOperator Ascaled( storage, minLevel, maxLevel );

   // Scale Laplace operator with viscosity
   Ascaled.scale( viscosity );

   hhg::P1DivxOperator          div_x( storage, minLevel, maxLevel );
   hhg::P1DivyOperator          div_y( storage, minLevel, maxLevel );
   hhg::P1DivTxOperator         divT_x( storage, minLevel, maxLevel );
   hhg::P1DivTyOperator         divT_y( storage, minLevel, maxLevel );
   hhg::P1LumpedInvMassOperator invDiagMass( storage, minLevel, maxLevel );

   std::array< std::shared_ptr< hhg::P1Function< real_t > >, 2 > velocity{
       {std::shared_ptr< hhg::P1Function< real_t > >( &u, boost::null_deleter() ),
        std::shared_ptr< hhg::P1Function< real_t > >( &v, boost::null_deleter() )}};
   hhg::DGUpwindOperator< hhg::P1Function< real_t > > N( storage, velocity, minLevel, maxLevel );

   typedef hhg::CGSolver< hhg::P1Function< real_t >, hhg::P1LaplaceOperator > CoarseSolver;
   typedef P1toP1LinearRestriction RestrictionOperator;
   typedef P1toP1LinearProlongation ProlongationOperator;

   auto coarseLaplaceSolver = std::make_shared< CoarseSolver >( storage, minLevel, minLevel );
   RestrictionOperator restrictionOperator;
   ProlongationOperator prolongationOperator;

   typedef GMultigridSolver< hhg::P1Function< real_t >, hhg::P1LaplaceOperator, CoarseSolver, RestrictionOperator, ProlongationOperator > LaplaceSover;
   LaplaceSover laplaceSolver( storage, coarseLaplaceSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   u.interpolate( bc_x, maxLevel, hhg::DirichletBoundary );
   v.interpolate( bc_y, maxLevel, hhg::DirichletBoundary );
   p.interpolate( zero, maxLevel - 1, hhg::NeumannBoundary );
   ones.interpolate( one, maxLevel, hhg::All );

   u_dg->projectP1( u, maxLevel, hhg::All );
   v_dg->projectP1( v, maxLevel, hhg::All );

   hhg::VTKOutput vtkOutput( "../output", "test", plotModulo );
   vtkOutput.add( &u );
   vtkOutput.add( &v );
   vtkOutput.add( &p );
   vtkOutput.write( maxLevel, iter );
   ++iter;

   while( time < endTime )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "time = " << time );
      time += dt;
      u.interpolate( bc_x, maxLevel, hhg::DirichletBoundary );

      u_dg_old->projectP1( u, maxLevel, hhg::All );
      v_dg_old->projectP1( v, maxLevel, hhg::All );

      N.apply( *u_dg_old, *u_dg, maxLevel, hhg::All, Replace );
      N.apply( *v_dg_old, *v_dg, maxLevel, hhg::All, Replace );

      // Predict u
      tmp.integrateDG( *u_dg, ones, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      Ascaled.apply( u, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary, Add );
      invDiagMass.apply( tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      u.add( {-dt}, {&tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary );

      // Predict v
      tmp.integrateDG( *v_dg, ones, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      Ascaled.apply( v, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary, Add );
      invDiagMass.apply( tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      v.add( {-dt}, {&tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary );

      // Solve p
      p.interpolate( zero, maxLevel - 1, hhg::NeumannBoundary );
      div_x.apply( u, p_rhs, maxLevel, hhg::Inner | hhg::DirichletBoundary, Replace );
      div_y.apply( v, p_rhs, maxLevel, hhg::Inner | hhg::DirichletBoundary, Add );

      restrictionOperator( p_rhs, maxLevel, hhg::Inner | hhg::DirichletBoundary );

      if( !neumann )
      {
         hhg::vertexdof::projectMean( p_rhs, tmp, maxLevel - 1 );
      }

      for( uint_t outer = 0; outer < outerIterations; ++outer )
      {
         laplaceSolver.solve( A,
                              p,
                              p_rhs,
                              p_res,
                              maxLevel - 1,
                              1e-2,
                              max_cg_iter,
                              hhg::Inner | hhg::DirichletBoundary,
                              LaplaceSover::CycleType::VCYCLE,
                              false );
      }

      if( !neumann )
      {
         hhg::vertexdof::projectMean( p, tmp, maxLevel - 1 );
      }

      prolongationOperator( p, maxLevel - 1, hhg::Inner | hhg::DirichletBoundary );

      // Correct u
      divT_x.apply( p, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      invDiagMass.apply( tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      u.add( {-1.0}, {&tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary );

      // Correct v
      divT_y.apply( p, tmp, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      invDiagMass.apply( tmp, tmp2, maxLevel, hhg::Inner | hhg::NeumannBoundary );
      v.add( {-1.0}, {&tmp2}, maxLevel, hhg::Inner | hhg::NeumannBoundary );

      vtkOutput.write( maxLevel, iter );
      ++iter;

      u_dg_old.swap( u_dg );
      v_dg_old.swap( v_dg );
   }

   timingTree->stop( "Global" );
   auto reduced_tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( reduced_tt );

   return 0;
}
