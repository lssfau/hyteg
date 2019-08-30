#include <cmath>

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"

using walberla::real_c;
using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = "../../data/meshes/quad_4el_neumann.msh";

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 2;
   //size_t maxiter  = 10000;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );

   hyteg::P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > bc_x = []( const hyteg::Point3D& x ) {
      if( x[0] < 1e-8 )
      {
         return 16.0 * ( x[1] - 0.5 ) * ( 1.0 - x[1] );
      } else
      {
         return 0.0;
      }
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.u.interpolate( bc_x, maxLevel, hyteg::DirichletBoundary );
   u.v.interpolate( zero, maxLevel, hyteg::DirichletBoundary );

   auto solver = hyteg::MinResSolver< hyteg::P1StokesOperator >( storage, minLevel, maxLevel );
   solver.solve( L, u, f, maxLevel );

   L.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
   real_t final_residuum = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) ) /
                           real_c( hyteg::numberOfGlobalDoFs< hyteg::P1StokesFunctionTag >( *storage, maxLevel ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residuum: " << final_residuum )

   WALBERLA_CHECK_LESS( final_residuum, 9.1e-07 );
   //hyteg::VTKWriter<hyteg::P1Function< real_t >>({ u.u, u.v, u.p }, maxLevel, "../output", "stokes_stab");
   return EXIT_SUCCESS;
}
