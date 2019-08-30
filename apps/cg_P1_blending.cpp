#include <core/Environment.h>
#include <core/config/Config.h>

#include "core/timing/Timer.h"

#include "hyteg/VTKWriter.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main( int argc, char* argv[] )
{
   /// create enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const size_t level   = 4;
   const uint_t maxiter = 10000;

   /// read mesh file and create storage
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../data/meshes/unitsquare_with_circular_hole.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{{0.5, 0.5, 0}};
   real_t  circleRadius = 0.25;

   for ( const auto& it : setupStorage.getFaces() )
   {
      Face& face = *( it.second );

      std::vector< PrimitiveID > neighborEdgesOnBoundary = face.neighborEdges();
      neighborEdgesOnBoundary.erase(
          std::remove_if( neighborEdgesOnBoundary.begin(),
                          neighborEdgesOnBoundary.end(),
                          [&setupStorage]( const PrimitiveID& id ) { return !setupStorage.onBoundary( id ); } ),
          neighborEdgesOnBoundary.end() );

      if ( neighborEdgesOnBoundary.size() > 0 )
      {
         Edge& edge = *setupStorage.getEdge( neighborEdgesOnBoundary[0] );

         if ( ( edge.getCoordinates()[0] - circleCenter ).norm() < 0.4 )
         {
            setupStorage.setGeometryMap( edge.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
            setupStorage.setGeometryMap( face.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
         }
      }
   }

   hhg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   typedef hhg::P1BlendingLaplaceOperator SolveOperator;
   SolveOperator                          L( storage, level, level );

   P1BlendingMassOperator M( storage, level, level );

   auto x       = std::make_shared< hhg::P1Function< real_t > >( "x", storage, level, level );
   auto y       = std::make_shared< hhg::P1Function< real_t > >( "y", storage, level, level );
   auto u       = std::make_shared< hhg::P1Function< real_t > >( "u", storage, level, level );
   auto u_exact = std::make_shared< hhg::P1Function< real_t > >( "u_exact", storage, level, level );
   auto f       = std::make_shared< hhg::P1Function< real_t > >( "f", storage, level, level );
   auto r       = std::make_shared< hhg::P1Function< real_t > >( "r", storage, level, level );
   auto err     = std::make_shared< hhg::P1Function< real_t > >( "err", storage, level, level );
   auto helper  = std::make_shared< hhg::P1Function< real_t > >( "helper", storage, level, level );

   std::function< real_t( const hhg::Point3D& ) > ones  = []( const hhg::Point3D& ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > zeros = []( const hhg::Point3D& ) { return 0.0; };

   helper->interpolate( ones, level );
   real_t npoints = helper->dotGlobal( *helper, level );

   std::function< real_t( const hhg::Point3D& ) > tmp_x = [&]( const hhg::Point3D& x_ ) { return x_[0]; };

   std::function< real_t( const hhg::Point3D& ) > tmp_y = [&]( const hhg::Point3D& x_ ) { return x_[1]; };

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x_ ) {
      return sin( x_[0] ) * cos( x_[1] ) / ( x_[0] * x_[1] + 1 );
   };
   std::function< real_t( const hhg::Point3D& ) > rhs = []( const hhg::Point3D& x_ ) {
      return 2 *
             ( -( pow( x_[0], 2 ) + pow( x_[1], 2 ) ) * sin( x_[0] ) * cos( x_[1] ) +
               pow( x_[0] * x_[1] + 1, 2 ) * sin( x_[0] ) * cos( x_[1] ) +
               ( x_[0] * x_[1] + 1 ) * ( -x_[0] * sin( x_[0] ) * sin( x_[1] ) + x_[1] * cos( x_[0] ) * cos( x_[1] ) ) ) /
             pow( x_[0] * x_[1] + 1, 3 );
   };

   x->interpolate( tmp_x, level, hhg::All );
   y->interpolate( tmp_y, level, hhg::All );

   u->interpolate( exact, level, hhg::DirichletBoundary );
   u_exact->interpolate( exact, level );
   helper->interpolate( rhs, level, hhg::All );
   M.apply( *helper, *f, level, hhg::All );

   auto solver = hhg::CGSolver< SolveOperator >( storage, level, level, maxiter, 1e-10 );
   solver.solve( L, *u, *f, level );

   err->assign( {1.0, -1.0}, {*u, *u_exact}, level, hhg::All );

   real_t discr_l2_err = std::sqrt( err->dotGlobal( *err, level ) / npoints );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   VTKOutput vtkOutput( "../output", "cg_P1_blending", storage );
   vtkOutput.add( *x );
   vtkOutput.add( *y );
   vtkOutput.add( *u.get() );
   vtkOutput.add( *u_exact );
   vtkOutput.add( *err );
   vtkOutput.write( level );

   return 0;
}
