#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"

#include "tinyhhg_core/likwidwrapper.hpp"

#include "core/Environment.h"

using walberla::real_t;
using walberla::real_c;
using namespace hhg;


class exactFunctor {
public:
  virtual real_t operator()(const hhg::Point3D& x) const = 0;
};


template< typename ValueType, uint_t Level >
inline void interpolateStdFunction(Face &face,
                            const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& faceMemoryId,
                            std::function<ValueType(const hhg::Point3D &)> &expr) {

  FaceP1FunctionMemory< ValueType > *faceMemory = face.getData(faceMemoryId);
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;
  auto dstPtr = faceMemory->getPointer( Level );
  x0 = face.coords[0];
  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += real_c(i)*d2 + d0;

    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      dstPtr[vertexdof::macroface::indexFromVertex<Level>(j, i, stencilDirection::VERTEX_C)] = expr(x);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

template< typename ValueType, uint_t Level, typename Expr >
inline void interpolateTemplate(Face &face,
                         const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& faceMemoryId,
                         const Expr& expr) {
  FaceP1FunctionMemory< ValueType > *faceMemory = face.getData(faceMemoryId);
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;
  auto dstPtr = faceMemory->getPointer( Level );
  x0 = face.coords[0];
  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += real_c(i)*d2 + d0;

    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      dstPtr[vertexdof::macroface::indexFromVertex<Level>(j, i, stencilDirection::VERTEX_C)] = expr(x);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

template< typename ValueType, uint_t Level >
inline void interpolateFunctor(Face &face,
                                const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& faceMemoryId,
                                const exactFunctor& exprFunctor) {

  FaceP1FunctionMemory< ValueType > *faceMemory = face.getData(faceMemoryId);
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;
  auto dstPtr = faceMemory->getPointer( Level );
  x0 = face.coords[0];
  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));
  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += real_c(i)*d2 + d0;

    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      dstPtr[vertexdof::macroface::indexFromVertex<Level>(j, i, stencilDirection::VERTEX_C)] = exprFunctor(x);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

template< typename ValueType, uint_t Level >
inline void interpolateWithoutFunction(Face &face,
                                       const PrimitiveDataID<FaceP1FunctionMemory< ValueType >, Face>& faceMemoryId) {
  FaceP1FunctionMemory< ValueType > *faceMemory = face.getData(faceMemoryId);
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x, x0;
  auto dstPtr = faceMemory->getPointer( Level );
  x0 = face.coords[0];
  Point3D d0 = (face.coords[1] - face.coords[0])/(walberla::real_c(rowsize - 1));
  Point3D d2 = (face.coords[2] - face.coords[0])/(walberla::real_c(rowsize - 1));

  uint_t inner_rowsize = rowsize;

  for (uint_t i = 1; i < rowsize - 2; ++i) {
    x = x0;
    x += real_c(i)*d2 + d0;

    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      dstPtr[vertexdof::macroface::indexFromVertex<Level>(j, i, stencilDirection::VERTEX_C)] = sqrt(x[0] * x[0] + x[1] * x[1]);
      x += d0;
    }

    inner_rowsize -= 1;
  }
}

class derivedFunctor : public exactFunctor{
public:
  real_t operator()(const hhg::Point3D& x) const{
    return sqrt(x[0] * x[0] + x[1] * x[1]);
  }
};

real_t exact(const hhg::Point3D& x) {
  return sqrt(x[0] * x[0] + x[1] * x[1]);
}

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

  auto x = std::make_shared<hhg::P1Function<real_t>>("x", storage, level, level);
  std::shared_ptr<Face> face = storage->getFaces().begin().operator*().second;


  std::function<real_t(const hhg::Point3D&)> exactFunc =
    [&](const hhg::Point3D& point) { return sqrt(point[0] * point[0] + point[1] * point[1]); };

  //P1Function< real_t > x("x", storage, level, level);

  walberla::WcTimer timer;

  LIKWID_MARKER_START("std::function");
  timer.reset();
  interpolateStdFunction< real_t, level >(*face,x->getFaceDataID(), exactFunc);
  timer.end();
  LIKWID_MARKER_STOP("std::function");
  WALBERLA_LOG_INFO_ON_ROOT( std::setw(20) << "std::function: " <<  timer.last() << " " << x->dot(*x, level, hhg::Inner));

  LIKWID_MARKER_START("Template");
  timer.reset();
  interpolateTemplate< real_t, level >(*face,x->getFaceDataID(), exact);
  timer.end();
  LIKWID_MARKER_STOP("Template");
  WALBERLA_LOG_INFO_ON_ROOT(std::setw(20) << "Template: " <<  timer.last() << " " << x->dot(*x, level, hhg::Inner));

  LIKWID_MARKER_START("without Function");
  timer.reset();
  interpolateWithoutFunction< real_t, level >(*face,x->getFaceDataID());
  timer.end();
  LIKWID_MARKER_STOP("without Function");
  WALBERLA_LOG_INFO_ON_ROOT(std::setw(20) << "Without Function: " << timer.last() << " " << x->dot(*x, level, hhg::Inner));

  LIKWID_MARKER_START("Functor");
  timer.reset();
  interpolateFunctor< real_t, level >(*face,x->getFaceDataID(),derivedFunctor());
  timer.end();
  LIKWID_MARKER_STOP("Functor");
  WALBERLA_LOG_INFO_ON_ROOT(std::setw(20) << "Functor: " << timer.last() << " " << x->dot(*x, level, hhg::Inner));



  LIKWID_MARKER_CLOSE;




}
