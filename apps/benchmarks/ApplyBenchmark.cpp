#include <tinyhhg_core/tinyhhg.hpp>
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"

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

  const size_t level = 8;

  auto src = std::make_shared<hhg::P1Function<real_t>>("src", storage, level, level);
  auto dst = std::make_shared<hhg::P1Function<real_t>>("dst", storage, level, level);

  hhg::P1MassOperator M(storage, level, level);


  std::shared_ptr<Face> face = storage->getFaces().begin().operator*().second;


  std::function<real_t(const hhg::Point3D&)> exactFunc =
    [&](const hhg::Point3D& point) { return sqrt(point[0] * point[0] + point[1] * point[1]); };

  //P1Function< real_t > x("x", storage, level, level);
  src->interpolate(exactFunc,level);

  walberla::WcTimer timer;

  LIKWID_MARKER_START("apply");
  timer.reset();
  vertexdof::macroface::apply_tmpl<real_t, level >(*face,M.getFaceStencilID(),src->getFaceDataID(),dst->getFaceDataID(),Replace);
  timer.end();
  LIKWID_MARKER_STOP("apply");
  WALBERLA_LOG_INFO_ON_ROOT("time with walberla timer: " << timer.last() );

  /// do something with the result to prevent the compiler from removing all the computations
  real_t check = vertexdof::macroface::dotTmpl< real_t, level >(*face,dst->getFaceDataID(),dst->getFaceDataID());
  WALBERLA_CHECK_FLOAT_UNEQUAL(check ,0. );

  LIKWID_MARKER_CLOSE;




}
