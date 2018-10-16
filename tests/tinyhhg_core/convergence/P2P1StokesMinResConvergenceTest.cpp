#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"

using walberla::real_t;
using walberla::real_c;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = "../../data/meshes/quad_4el_neumann.msh";

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   size_t level   = 2;
   size_t maxiter = 10000;

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::P2P1TaylorHoodFunction< real_t > r( "r", storage, level, level );
   hhg::P2P1TaylorHoodFunction< real_t > f( "f", storage, level, level );
   hhg::P2P1TaylorHoodFunction< real_t > u( "u", storage, level, level );

   hhg::P2P1TaylorHoodStokesOperator L( storage, level, level );

   std::function< real_t( const hhg::Point3D& ) > bc_x = []( const hhg::Point3D& x ) { return 4.0 * ( 1.0 - x[1] ) * x[1]; };
   std::function< real_t( const hhg::Point3D& ) > rhs  = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   u.u.interpolate( zero, level );
   u.u.interpolate( bc_x, level, hhg::DirichletBoundary );
   u.v.interpolate( zero, level, hhg::DirichletBoundary );

   auto solver =
       hhg::MinResSolver< hhg::P2P1TaylorHoodFunction< real_t >, hhg::P2P1TaylorHoodStokesOperator >( storage, level, level );
   solver.solve( L, u, f, r, level, 1e-3, maxiter, hhg::Inner | hhg::NeumannBoundary, true );

   L.apply( u, r, level, hhg::Inner | hhg::NeumannBoundary );
   real_t final_residuum = std::sqrt( r.dotGlobal( r, level, hhg::Inner ) ) /
                           real_c( hhg::numberOfGlobalDoFs< hhg::P1StokesFunctionTag >( *storage, level ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Residuum: " << final_residuum )

   WALBERLA_CHECK_LESS( final_residuum, 5e-05 );
   return EXIT_SUCCESS;
}
