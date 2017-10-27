#include <tinyhhg_core/tinyhhg.hpp>
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1FaceIndex.hpp"

#include "tinyhhg_core/likwidwrapper.hpp"
#include "tinyhhg_core/indexing/Optimization.hpp"
#include "tinyhhg_core/indexing/VertexDoFIndexing.hpp"

#include "core/Environment.h"

using walberla::real_t;
using walberla::real_c;
using namespace hhg;

std::shared_ptr<PrimitiveStorage> globalStorage;

template<size_t Level>
inline void manualApply(real_t * oprPtr,
                        real_t* srcPtr,
                        real_t* dstPtr,
                        UpdateType update) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  size_t inner_rowsize = rowsize;

  size_t br = 1;
  size_t mr = 1 + rowsize;
  size_t tr = mr + (rowsize - 1);

  for (size_t i = 0; i < rowsize - 3; ++i) {
    for (size_t j = 0; j < inner_rowsize - 3; ++j) {
      if (update == Replace) {
        dstPtr[mr] = oprPtr[0] * srcPtr[br] + oprPtr[1] * srcPtr[br + 1]
                  + oprPtr[2] * srcPtr[mr - 1] + oprPtr[3] * srcPtr[mr] + oprPtr[4] * srcPtr[mr + 1]
                  + oprPtr[5] * srcPtr[tr - 1] + oprPtr[6] * srcPtr[tr];
      } else if (update == Add) {
        dstPtr[mr] += oprPtr[0] * srcPtr[br] + oprPtr[1] * srcPtr[br + 1]
                   + oprPtr[2] * srcPtr[mr - 1] + oprPtr[3] * srcPtr[mr] + oprPtr[4] * srcPtr[mr + 1]
                   + oprPtr[5] * srcPtr[tr - 1] + oprPtr[6] * srcPtr[tr];
      }

      br += 1;
      mr += 1;
      tr += 1;
    }

    br += 3;
    mr += 2;
    tr += 1;
    --inner_rowsize;
  }
}
template< uint_t Level >
inline void apply(real_t * oprPtr,
                        real_t* srcPtr,
                        real_t* dstPtr,
                       UpdateType update) {
  using namespace hhg::P1Face::FaceCoordsVertex;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  real_t tmp;

#ifdef WALBERLA_CXX_COMPILER_IS_INTEL
#pragma ivdep
#endif
  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = oprPtr[VERTEX_C]*srcPtr[index<Level>(i, j, VERTEX_C)];

      tmp += oprPtr[neighbors[0]]*srcPtr[index<Level>(i, j, neighbors[1])];
      tmp += oprPtr[neighbors[1]]*srcPtr[index<Level>(i, j, neighbors[1])];
      tmp += oprPtr[neighbors[2]]*srcPtr[index<Level>(i, j, neighbors[2])];
      tmp += oprPtr[neighbors[3]]*srcPtr[index<Level>(i, j, neighbors[3])];
      tmp += oprPtr[neighbors[4]]*srcPtr[index<Level>(i, j, neighbors[4])];
      tmp += oprPtr[neighbors[5]]*srcPtr[index<Level>(i, j, neighbors[5])];
      //for (auto neighbor : neighbors) {
//      for(uint_t k = 0; k < neighbors.size(); ++k){
//        tmp += oprPtr[neighbors[k]]*srcPtr[index<Level>(i, j, neighbors[k])];
//      }

      if (update==Replace) {
        dstPtr[index<Level>(i, j, VERTEX_C)] = tmp;
      } else if (update==Add) {
        dstPtr[index<Level>(i, j, VERTEX_C)] += tmp;
      }
    }
    --inner_rowsize;
  }
}

template< uint_t Level >
inline void applyOptimized(real_t * oprPtr,
                  real_t* srcPtr,
                  real_t* dstPtr,
                  UpdateType update) {
  using namespace hhg::P1Face::FaceCoordsVertex;
  using namespace hhg::indexing;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  real_t tmp;

  uint_t actualCol;
  uint_t actualRow;
#ifdef WALBERLA_CXX_COMPILER_IS_INTEL
#pragma ivdep
#endif
  for ( const auto & it : FaceIterator< vertexdof::levelToWidth< Level >, 1 >() ){

    actualCol = it.col();
    actualRow = it.row();

    tmp = oprPtr[VERTEX_C]*srcPtr[index<Level>(actualCol, actualRow, VERTEX_C)];

    tmp += oprPtr[neighbors[0]]*srcPtr[index<Level>(actualCol, actualRow, neighbors[1])];
    tmp += oprPtr[neighbors[1]]*srcPtr[index<Level>(actualCol, actualRow, neighbors[1])];
    tmp += oprPtr[neighbors[2]]*srcPtr[index<Level>(actualCol, actualRow, neighbors[2])];
    tmp += oprPtr[neighbors[3]]*srcPtr[index<Level>(actualCol, actualRow, neighbors[3])];
    tmp += oprPtr[neighbors[4]]*srcPtr[index<Level>(actualCol, actualRow, neighbors[4])];
    tmp += oprPtr[neighbors[5]]*srcPtr[index<Level>(actualCol, actualRow, neighbors[5])];

    if (update==Replace) {
      dstPtr[index<Level>(actualCol, actualRow, VERTEX_C)] = tmp;
    } else if (update==Add) {
      dstPtr[index<Level>(actualCol, actualRow, VERTEX_C)] += tmp;
    }
  }

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
  globalStorage = storage;

  const size_t level = 14;

  auto src = std::make_shared<hhg::P1Function<real_t>>("src", storage, level, level);
  auto dst1 = std::make_shared<hhg::P1Function<real_t>>("dst", storage, level, level);
  auto dst2 = std::make_shared<hhg::P1Function<real_t>>("dst", storage, level, level);


  hhg::P1LaplaceOperator M(storage, level, level);

  std::shared_ptr<Face> face = storage->getFaces().begin().operator*().second;
  //auto &oprPtr2 = face->getData(M.getFaceStencilID())->data[level];

  std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D& x)  { return x.x[0] * 4; };
  src->interpolate(ones,level);

  real_t* oprPtr = face->getData(M.getFaceStencilID())->data[level].get();
  for(uint_t i = 0; i < 7; ++i){
    oprPtr[i] = real_c(i);
  }
  real_t* srcPtr = face->getData(src->getFaceDataID())->getPointer( level );
  real_t* dst1Ptr = face->getData(dst1->getFaceDataID())->getPointer( level );
  real_t* dst2Ptr = face->getData(dst2->getFaceDataID())->getPointer( level );

  walberla::WcTimer timer;
  timer.reset();
  LIKWID_MARKER_START("Optimized");
  applyOptimized< level >(oprPtr,srcPtr,dst1Ptr,Replace);
  LIKWID_MARKER_STOP("Optimized");
  timer.end();
  WALBERLA_LOG_INFO_ON_ROOT(std::setw(20) << "Optimized: " << timer.last());

  timer.reset();
  LIKWID_MARKER_START("Index");
  apply< level >(oprPtr,srcPtr,dst2Ptr,Replace);
  LIKWID_MARKER_STOP("Index");
  timer.end();
  WALBERLA_LOG_INFO_ON_ROOT(std::setw(20) << "Index: " << timer.last());


  real_t sum1 = 0,sum2 = 0;

  for(uint_t i = 0; i < levelinfo::num_microvertices_per_face(level); ++i){
    WALBERLA_CHECK_FLOAT_EQUAL(dst1Ptr[i],dst2Ptr[i], "i was: " << i );
    sum1 += dst1Ptr[i];
    sum2 += dst2Ptr[i];
  }

  WALBERLA_LOG_INFO_ON_ROOT("sum1: " << sum1 << " sum2: " << sum2);

  LIKWID_MARKER_CLOSE;

}