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

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 2;
   size_t maxiter  = 10000;

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::P1StokesFunction< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );

   hhg::P1StokesOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > bc_x = []( const hhg::Point3D& x ) {
      if( x[0] < 1e-8 )
      {
         return 16.0 * ( x[1] - 0.5 ) * ( 1.0 - x[1] );
      } else
      {
         return 0.0;
      }
   };
   std::function< real_t( const hhg::Point3D& ) > rhs  = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   u.u.interpolate( bc_x, maxLevel, hhg::DirichletBoundary );
   u.v.interpolate( zero, maxLevel, hhg::DirichletBoundary );

   auto solver = hhg::MinResSolver< hhg::P1StokesFunction< real_t >, hhg::P1StokesOperator >( storage, minLevel, maxLevel );
   solver.solve( L, u, f, r, maxLevel, 1e-5, maxiter, hhg::Inner | hhg::NeumannBoundary, true );

   L.apply( u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary );
   real_t final_residuum = std::sqrt( r.dot( r, maxLevel, hhg::Inner ) ) /
                           real_c( hhg::numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, maxLevel ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residuum: " << final_residuum )

   WALBERLA_CHECK_LESS( final_residuum, 8.1e-07 );
   //hhg::VTKWriter<hhg::P1Function< real_t >>({ &u.u, &u.v, &u.p }, maxLevel, "../output", "stokes_stab");
   return EXIT_SUCCESS;
}
