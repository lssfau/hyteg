
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"
#include "core/debug/all.h"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using namespace hhg;

int main(int argc, char **argv)
{
  walberla::debug::enterTestMode();
  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/3D/tet_1el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  const uint_t level = 2;

  EdgeDoFFunction< real_t > x("x", storage, level, level);

  std::function< real_t( const Point3D & )> exact = []( const Point3D & xx )
  {
    return 2*xx[0] + xx[1] * xx[2];
  };

  VTKOutput vtkOutput( "../../output", "EdgeDoFInterpolation3DTest" );
  vtkOutput.set3D();
  vtkOutput.add( &x );

  vtkOutput.write( level, 0 );

  x.interpolate( exact, level );

  vtkOutput.write( level, 1 );
}
