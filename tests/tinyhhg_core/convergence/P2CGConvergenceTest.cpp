#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hhg {

void P2CGTest(const std::string &meshFile, const uint_t level, const real_t targetError, const bool localMPI)
{
   const auto meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2CGConvergenceTest_domain" );

   hhg::P2ConstantLaplaceOperator L( storage, level, level );

   hhg::P2Function< real_t > r( "r", storage, level, level );
   hhg::P2Function< real_t > f( "f", storage, level, level );
   hhg::P2Function< real_t > u( "u", storage, level, level );
   hhg::P2Function< real_t > u_exact( "u_exact", storage, level, level );
   hhg::P2Function< real_t > err( "err", storage, level, level );
   hhg::P2Function< real_t > npoints_helper( "npoints_helper", storage, level, level );

   if(localMPI){
      u.setLocalCommunicationMode( communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
   }

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hhg::Point3D& ) > rhs   = []( const hhg::Point3D& ) { return 0; };
   std::function< real_t( const hhg::Point3D& ) > ones  = []( const hhg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hhg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   auto solver = hhg::CGSolver< hhg::P2ConstantLaplaceOperator >( storage, level, level );
   solver.solve( L, u, f, level );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );
   npoints_helper.interpolate( ones, level );

   const real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   hhg::VTKOutput vtkOutput( "../../output", "P2CGConvergenceTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( f );
   vtkOutput.add( r );
   vtkOutput.add( err );
   vtkOutput.add( npoints_helper );
   vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_CHECK_LESS( discr_l2_err, targetError );
}

} // namespace hhg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hhg::P2CGTest("../../data/meshes//tri_1el.msh", 3, 1e-7, false);
   hhg::P2CGTest("../../data/meshes//quad_4el.msh", 3, 1e-6, false);
   hhg::P2CGTest("../../data/meshes/3D/tet_1el.msh", 2, 3e-6, false);
   hhg::P2CGTest("../../data/meshes/3D/tet_1el.msh", 3, 3e-7, true);
   hhg::P2CGTest("../../data/meshes/3D/pyramid_2el.msh", 2, 3e-5, false);
   hhg::P2CGTest("../../data/meshes/3D/regular_octahedron_8el.msh", 2, 1.7e-5, true);
}