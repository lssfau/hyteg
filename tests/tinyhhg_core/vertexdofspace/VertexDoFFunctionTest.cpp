
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/vtkwriter.hpp"

namespace hhg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;

static void testVertexDoFFunction()
{
  const uint_t level = 3;

  MeshInfo mesh  = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  writeDomainPartitioningVTK( storage, "../../output", "single_tet" );

  auto x = std::make_shared< P1Function< real_t > >( "x", storage, level, level );

  VTKOutput vtkOutput( "../../output", "vertex_dof_macro_cell_test" );
  vtkOutput.set3D();
  vtkOutput.add( x );
  vtkOutput.write( level );

}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testVertexDoFFunction();

   return EXIT_SUCCESS;
}
