#include "core/DataTypes.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P1BlendingStokesOperator.hpp"
#include "tinyhhg_core/composites/P1CoefficientStokesOperator.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/geometry/CircularMap.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/unitsquare_with_circular_hole_neumann.msh";

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{{0.5, 0.5, 0}};
   real_t  circleRadius = 0.25;

   for( auto it = setupStorage.beginFaces(); it != setupStorage.endFaces(); ++it )
   {
      Face& face = *it->second;

      if( face.hasBoundaryEdge() )
      {
         Edge& edge = *setupStorage.getEdge( face.edgesOnBoundary[0] );

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

   size_t minLevel = 2;
   size_t maxLevel = 4;
   size_t maxiter  = 10000;

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::P1Function< real_t >                    tmp( "tmp", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              u( "u", storage, minLevel, maxLevel );
   std::shared_ptr< hhg::P1Function< real_t > > coefficient =
       std::make_shared< hhg::P1Function< real_t > >( "coeff", storage, minLevel, maxLevel );

   typedef hhg::P1BlendingStokesOperator SolveOperator;

   SolveOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > coeff  = []( const hhg::Point3D& x ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > zero   = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones   = []( const hhg::Point3D& ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > inflow = []( const hhg::Point3D& x ) {
      if( x[0] < 1e-4 )
      {
         return 4.0 * x[1] * ( 1.0 - x[1] );
      }

      return 0.0;
   };

   u.u.interpolate( inflow, maxLevel, hhg::DirichletBoundary );
   u.v.interpolate( zero, maxLevel, hhg::DirichletBoundary );
   coefficient->interpolate( coeff, maxLevel );

   auto solver = hhg::MinResSolver< hhg::P1StokesFunction< real_t >, SolveOperator >( storage, minLevel, maxLevel );
   solver.solve( L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner | hhg::NeumannBoundary, true );

   // u_u*iHat + u_v*jHat
   hhg::VTKOutput vtkOutput( "../output", "stokes_stab_varcoeff" );
   vtkOutput.add( &u.u );
   vtkOutput.add( &u.v );
   vtkOutput.add( &u.p );
   vtkOutput.write( maxLevel, 0 );
   return EXIT_SUCCESS;
}
