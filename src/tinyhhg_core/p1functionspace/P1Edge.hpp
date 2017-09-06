#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "P1Memory.hpp"
#include "P1EdgeIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

namespace hhg {

namespace P1Edge {

template<uint_t Level>
inline void interpolateTmpl(Edge &edge,
                        const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &edgeMemoryId,
                        std::function<real_t(const hhg::Point3D &)> &expr) {
  using namespace EdgeCoordsVertex;

  EdgeP1FunctionMemory *edgeMemory = edge.getData(edgeMemoryId);

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x = edge.getCoordinates()[0];
  Point3D dx = edge.getDirection()/(real_t) (rowsize - 1);

  x += dx;

  for (size_t i = 1; i < rowsize - 1; ++i) {
    edgeMemory->data[Level][index<Level>(i, VERTEX_C)] = expr(x);
    x += dx;
  }
}

SPECIALIZE(void, interpolateTmpl, interpolate)

template<uint_t Level>
inline void assignTmpl(Edge &edge,
                   const std::vector<real_t> &scalars,
                   const std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>> &srcIds,
                   const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    real_t tmp = scalars[0]*edge.getData(srcIds[0])->data[Level][index<Level>(i, VERTEX_C)];

    for (size_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->data[Level][index<Level>(i, VERTEX_C)];
    }

    edge.getData(dstId)->data[Level][index<Level>(i, VERTEX_C)] = tmp;
  }
}

SPECIALIZE(void, assignTmpl, assign)

template<uint_t Level>
inline void addTmpl(Edge &edge,
                const std::vector<real_t> &scalars,
                const std::vector<PrimitiveDataID<EdgeP1FunctionMemory, Edge>> &srcIds,
                const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    real_t tmp = 0.0;

    for (size_t k = 0; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->data[Level][index<Level>(i, VERTEX_C)];
    }

    edge.getData(dstId)->data[Level][index<Level>(i, VERTEX_C)] += tmp;
  }
}

SPECIALIZE(void, addTmpl, add)

template<uint_t Level>
inline real_t dotTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &lhsMemoryId,
                  const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &rhsMemoryId) {
  using namespace EdgeCoordsVertex;

  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    sp += edge.getData(lhsMemoryId)->data[Level][index<Level>(i, VERTEX_C)]
        * edge.getData(rhsMemoryId)->data[Level][index<Level>(i, VERTEX_C)];
  }

  return sp;
}

SPECIALIZE(real_t, dotTmpl, dot)

template<uint_t Level>
inline void applyTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                  const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &srcId,
                  const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId, UpdateType update) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto &src = edge.getData(srcId)->data[Level];
  auto &dst = edge.getData(dstId)->data[Level];

  real_t tmp;

  for (size_t i = 1; i < rowsize - 1; ++i) {

    tmp = opr_data[VERTEX_C] * src[index<Level>(i, VERTEX_C)];

    for (auto& neighbor : neighbors_on_edge) {
      tmp += opr_data[neighbor] * src[index<Level>(i, neighbor)];
    }

    for (auto& neighbor : neighbors_south) {
      tmp += opr_data[neighbor] * src[index<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto& neighbor : neighbors_north) {
        tmp += opr_data[neighbor] * src[index<Level>(i, neighbor)];
      }
    }

    if (update == Replace) {
      dst[index<Level>(i, VERTEX_C)] = tmp;
    } else if (update == Add) {
      dst[index<Level>(i, VERTEX_C)] += tmp;
    }
  }
}

SPECIALIZE(void, applyTmpl, apply)

template<uint_t Level>
inline void smoothGSTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                      const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId,
                      const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &rhsId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto &dst = edge.getData(dstId)->data[Level];
  auto &rhs = edge.getData(rhsId)->data[Level];

  for (size_t i = 1; i < rowsize - 1; ++i) {

    dst[index<Level>(i, VERTEX_C)] = rhs[index<Level>(i, VERTEX_C)];

    for (auto& neighbor : neighbors_on_edge) {
      dst[index<Level>(i, VERTEX_C)] -= opr_data[neighbor] * dst[index<Level>(i, neighbor)];
    }

    for (auto& neighbor : neighbors_south) {
      dst[index<Level>(i, VERTEX_C)] -= opr_data[neighbor] * dst[index<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto& neighbor : neighbors_north) {
        dst[index<Level>(i, VERTEX_C)] -= opr_data[neighbor] * dst[index<Level>(i, neighbor)];
      }
    }

    dst[index<Level>(i, VERTEX_C)] /= opr_data[VERTEX_C];
  }
}

SPECIALIZE(void, smoothGSTmpl, smooth_gs)

template<uint_t Level>
inline void smoothJacTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                          const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId,
                          const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &rhsId,
                          const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &tmpId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto &dst = edge.getData(dstId)->data[Level];
  auto &rhs = edge.getData(rhsId)->data[Level];
  auto &tmp = edge.getData(tmpId)->data[Level];

  for (size_t i = 1; i < rowsize - 1; ++i) {

    dst[index<Level>(i, VERTEX_C)] = rhs[index<Level>(i, VERTEX_C)];

    for (auto& neighbor : neighbors_on_edge) {
      dst[index<Level>(i, VERTEX_C)] -= opr_data[neighbor] * tmp[index<Level>(i, neighbor)];
    }

    for (auto& neighbor : neighbors_south) {
      dst[index<Level>(i, VERTEX_C)] -= opr_data[neighbor] * tmp[index<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto& neighbor : neighbors_north) {
        dst[index<Level>(i, VERTEX_C)] -= opr_data[neighbor] * tmp[index<Level>(i, neighbor)];
      }
    }

    dst[index<Level>(i, VERTEX_C)] /= opr_data[VERTEX_C];
  }
}

SPECIALIZE(void, smoothJacTmpl, smooth_jac)

template<uint_t SourceLevel>
inline void prolongateTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &memoryId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize_c = levelinfo::num_microvertices_per_edge(SourceLevel);

  auto &edge_data_f = edge.getData(memoryId)->data[SourceLevel + 1];
  auto &edge_data_c = edge.getData(memoryId)->data[SourceLevel];

  size_t i_c;
  for (i_c = 1; i_c < rowsize_c - 1; ++i_c) {

    edge_data_f[index<SourceLevel+1>(2*i_c, VERTEX_C)] = edge_data_c[index<SourceLevel>(i_c, VERTEX_C)];
    edge_data_f[index<SourceLevel+1>(2*i_c-1, VERTEX_C)] = 0.5*(edge_data_c[index<SourceLevel>(i_c-1, VERTEX_C)]
                                                              + edge_data_c[index<SourceLevel>(i_c, VERTEX_C)]);
  }

  edge_data_f[index<SourceLevel+1>(2*i_c-1, VERTEX_C)] = 0.5*(edge_data_c[index<SourceLevel>(i_c-1, VERTEX_C)]
                                                            + edge_data_c[index<SourceLevel>(i_c, VERTEX_C)]);
}

SPECIALIZE(void, prolongateTmpl, prolongate)

template<uint_t Level>
inline void prolongateQuadraticTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &memoryId) {

  //TODO: rewrite using index function possible? maybe more generalized notion of Operator between different levels
  // is required.

  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(Level);
  size_t i_fine = 1;
  real_t invtemp = 1/8.;
  const real_t s1[3] = {3*invtemp, 6*invtemp, -invtemp};
  const real_t s2[3] = {-invtemp, 6*invtemp, 3*invtemp};

  auto &edge_data_f = edge.getData(memoryId)->data[Level + 1];
  auto &edge_data_c = edge.getData(memoryId)->data[Level];
  size_t i_coarse;
  for (i_coarse = 0; i_coarse < rowsize_coarse - 2; ++i_coarse) {
    edge_data_f[i_fine] =
        (s1[0]*edge_data_c[i_coarse] + s1[1]*edge_data_c[i_coarse + 1] + s1[2]*edge_data_c[i_coarse + 2]);
    edge_data_f[i_fine + 1] = edge_data_c[i_coarse + 1];
    i_fine += 2;
  }
  i_coarse--;
  edge_data_f[i_fine] =
      (s2[0]*edge_data_c[i_coarse] + s2[1]*edge_data_c[i_coarse + 1] + s2[2]*edge_data_c[i_coarse + 2]);

}

SPECIALIZE(void, prolongateQuadraticTmpl, prolongateQuadratic)

template<uint_t Level>
inline void restrictTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &memoryId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize_c = levelinfo::num_microvertices_per_edge(Level - 1);

  auto &edge_data_f = edge.getData(memoryId)->data[Level];
  auto &edge_data_c = edge.getData(memoryId)->data[Level - 1];

  uint_t i_c;
  for ( i_c = 1; i_c < rowsize_c - 1; ++i_c) {
    edge_data_c[index<Level-1>(i_c, VERTEX_C)] = 1.0 * edge_data_f[index<Level>(2*i_c, VERTEX_C)];

    for (auto& neighbor : neighbors_on_edge) {
      edge_data_c[index<Level-1>(i_c, VERTEX_C)] += 0.5 * edge_data_f[index<Level>(2*i_c, neighbor)];
    }

    for (auto& neighbor : neighbors_south) {
      edge_data_c[index<Level-1>(i_c, VERTEX_C)] += 0.5 * edge_data_f[index<Level>(2*i_c, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto& neighbor : neighbors_north) {
        edge_data_c[index<Level-1>(i_c, VERTEX_C)] += 0.5 * edge_data_f[index<Level>(2*i_c, neighbor)];
      }
    }
  }
}

SPECIALIZE(void, restrictTmpl, restrict)

template<uint_t Level>
inline void enumerateTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId, uint_t& num) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    edge.getData(dstId)->data[Level][index<Level>(i, VERTEX_C)] = walberla::real_c(num++);
  }
}

SPECIALIZE(void, enumerateTmpl, enumerate)

#ifdef HHG_BUILD_WITH_PETSC
template<uint_t Level>
inline void saveOperatorTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                         const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &srcId,
                         const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &dstId, Mat& mat) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto &src = edge.getData(srcId)->data[Level];
  auto &dst = edge.getData(dstId)->data[Level];


  for (uint_t i = 1; i < rowsize - 1; ++i) {
    PetscInt dstint = (PetscInt) dst[index<Level>(i, VERTEX_C)];
    PetscInt srcint = (PetscInt) src[index<Level>(i, VERTEX_C)];
    //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, VERTEX_C)], opr_data[VERTEX_C]);
    MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[VERTEX_C] ,INSERT_VALUES);         //TODO: Make this more efficient by grouping all of them in an array

    for (auto& neighbor : neighbors_on_edge) {
      srcint = (PetscInt) src[index<Level>(i, neighbor)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[neighbor] ,INSERT_VALUES);
    }

    for (auto& neighbor : neighbors_south) {
      srcint = (PetscInt) src[index<Level>(i, neighbor)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[neighbor] ,INSERT_VALUES);
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto& neighbor : neighbors_north) {
        srcint = (PetscInt) src[index<Level>(i, neighbor)];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[neighbor] ,INSERT_VALUES);
      }
    }
  }
}

SPECIALIZE(void, saveOperatorTmpl, saveOperator)

template<uint_t Level>
inline void createVectorFromFunctionTmpl(Edge &edge,
                                     const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &srcId,
                                     const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &numeratorId,
                                     Vec& vec) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &src = edge.getData(srcId)->data[Level];
  auto &numerator = edge.getData(numeratorId)->data[Level];

  PetscInt* numeratorInt = new PetscInt[rowsize-2];
  for(uint_t i = 0;i<rowsize-2; i++)
  {
    numeratorInt[i] = (PetscInt)numerator[i+1];
  }

  VecSetValues(vec,rowsize-2,numeratorInt,&src[1],INSERT_VALUES);
  delete[] numeratorInt;

}
SPECIALIZE(void, createVectorFromFunctionTmpl, createVectorFromFunction)

template<uint_t Level>
inline void createFunctionFromVectorTmpl(Edge &edge,
                                         const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &srcId,
                                         const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &numeratorId,
                                         Vec& vec) {
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &numerator = edge.getData(numeratorId)->data[Level];

  PetscInt* numeratorInt = new PetscInt[rowsize-2];
  for(uint_t i = 0;i<rowsize-2; i++)
  {
    numeratorInt[i] = (PetscInt)numerator[i+1];
  }

  VecGetValues(vec,rowsize-2,numeratorInt,&edge.getData(srcId)->data[Level][1]);
  delete[] numeratorInt;

}
SPECIALIZE(void, createFunctionFromVectorTmpl, createFunctionFromVector)

template<uint_t Level>
inline void applyDirichletBCTmpl(Edge &edge,std::vector<PetscInt> &mat,
                                 const PrimitiveDataID<EdgeP1FunctionMemory, Edge> &numeratorId){

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for(uint_t i = 1;i<rowsize-1; i++)
  {
    mat.push_back((PetscInt)edge.getData(numeratorId)->data[Level][i]);
  }

}
SPECIALIZE(void, applyDirichletBCTmpl, applyDirichletBC)
#endif

}
}
