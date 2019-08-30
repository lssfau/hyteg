#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/preconditioners/JacobiPreconditioner.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = "../../data/meshes/quad_4el.msh";

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 5;
   size_t maxiter  = 1000;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   hyteg::P1ConstantLaplaceOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) -> real_t {
      return x[0] * x[0] - x[1] * x[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs  = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, maxLevel, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );

   auto prec   = std::make_shared< hyteg::JacobiPreconditioner< hyteg::P1ConstantLaplaceOperator >>( storage, minLevel, maxLevel, 10 );
   auto solver = hyteg::MinResSolver< hyteg::P1ConstantLaplaceOperator >(
       storage, minLevel, maxLevel,maxiter,1e-8, prec );
   solver.solve( L, u, f, maxLevel );

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );

   npoints_helper.interpolate( ones, maxLevel );
   real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

   real_t discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << std::scientific << discr_l2_err );

   WALBERLA_CHECK_LESS( discr_l2_err, 6e-09 )

   //hyteg::VTKWriter<hyteg::P1Function< real_t >>({ u, u_exact, &f, &r, &err }, maxLevel, "../output", "minres");
   return EXIT_SUCCESS;
}
