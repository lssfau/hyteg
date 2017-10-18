#include "tinyhhg_core/dgfunctionspace/DGFunction.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

#include "tinyhhg_core/likwidwrapper.hpp"

#include "core/Environment.h"

using walberla::real_t;
using walberla::real_c;
using namespace hhg;

int main(int argc, char **argv) {

  LIKWID_MARKER_INIT;

  walberla::debug::enterTestMode();
  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::MPIManager::instance()->useWorldComm();

  LIKWID_MARKER_THREADINIT;

  MeshInfo meshInfo = MeshInfo::fromGmshFile("../data/meshes/tri_1el.msh");
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  const size_t level = 14;

  std::minstd_rand0 generator(42);

  std::function<real_t(const hhg::Point3D&)> exact = [&](const hhg::Point3D&) { return generator(); };

  std::shared_ptr<Face> face = storage->getFaces().begin().operator*().second;

  DGFunction< real_t > x("x", storage, level, level);

  LIKWID_MARKER_START("interpolate");
  DGFace::interpolate< real_t >(level,*face,x.getFaceDataID(),exact);
  LIKWID_MARKER_STOP("interpolate");




}