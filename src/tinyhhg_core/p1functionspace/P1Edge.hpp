#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/p1functionspace/P1EdgeIndex.hpp"
#include "tinyhhg_core/dgfunctionspace/DGEdgeIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

#include "core/DataTypes.h"

namespace hhg {

namespace P1Edge {

template<typename ValueType, uint_t Level>
inline ValueType assembleLocal(uint_t pos, const Matrix3r& localMatrix,
                               double* src,
                               double* coeff,
                               const std::array<EdgeCoordsVertex::DirVertex,3>& vertices,
                               const std::array<uint_t,3>& idx)
{
  using namespace EdgeCoordsVertex;

  ValueType meanCoeff = 1.0/3.0 * (coeff[index<Level>(pos, vertices[0])]
                                 + coeff[index<Level>(pos, vertices[1])]
                                 + coeff[index<Level>(pos, vertices[2])]);

  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[index<Level>(pos, vertices[0])]
         + localMatrix(idx[0],idx[1]) * src[index<Level>(pos, vertices[1])]
         + localMatrix(idx[0],idx[2]) * src[index<Level>(pos, vertices[2])];
  return meanCoeff * tmp;
}

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Edge &edge,
                        const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &edgeMemoryId,
                        std::function<ValueType(const hhg::Point3D &)> &expr) {
  using namespace EdgeCoordsVertex;

  EdgeP1FunctionMemory< ValueType > *edgeMemory = edge.getData(edgeMemoryId);

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  Point3D x = edge.getCoordinates()[0];
  Point3D dx = edge.getDirection()/(real_t) (rowsize - 1);

  x += dx;

  for (size_t i = 1; i < rowsize - 1; ++i) {
    edgeMemory->getPointer( Level )[index<Level>(i, VERTEX_C)] = expr(x);
    x += dx;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate )

template< typename ValueType, uint_t Level >
inline void assignTmpl(Edge &edge,
                   const std::vector<ValueType> &scalars,
                   const std::vector<PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge>> &srcIds,
                   const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    ValueType tmp = scalars[0]*edge.getData(srcIds[0])->getPointer( Level )[index<Level>(i, VERTEX_C)];

    for (size_t k = 1; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( Level )[index<Level>(i, VERTEX_C)];
    }

    edge.getData(dstId)->getPointer( Level )[index<Level>(i, VERTEX_C)] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, assignTmpl, assign )

template< typename ValueType, uint_t Level >
inline void addTmpl(Edge &edge,
                const std::vector<ValueType> &scalars,
                const std::vector<PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge>> &srcIds,
                const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    ValueType tmp = 0.0;

    for (size_t k = 0; k < srcIds.size(); ++k) {
      tmp += scalars[k]*edge.getData(srcIds[k])->getPointer( Level )[index<Level>(i, VERTEX_C)];
    }

    edge.getData(dstId)->getPointer( Level )[index<Level>(i, VERTEX_C)] += tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add )

template< typename ValueType, uint_t Level >
inline real_t dotTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &lhsMemoryId,
                  const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &rhsMemoryId) {
  using namespace EdgeCoordsVertex;

  real_t sp = 0.0;
  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    sp += edge.getData(lhsMemoryId)->getPointer( Level )[index<Level>(i, VERTEX_C)]
        * edge.getData(rhsMemoryId)->getPointer( Level )[index<Level>(i, VERTEX_C)];
  }

  return sp;
}

SPECIALIZE_WITH_VALUETYPE( real_t, dotTmpl, dot )

template< typename ValueType, uint_t Level >
inline void applyTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                  const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &srcId,
                  const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId, UpdateType update) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  ValueType tmp;

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

SPECIALIZE_WITH_VALUETYPE( void, applyTmpl, apply )

template< typename ValueType, uint_t Level >
inline void applyCoefficientTmpl(Edge &edge,
                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                 const PrimitiveDataID<EdgeP1LocalMatrixMemory, Edge> &operatorId,
                                 const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &srcId,
                                 const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId,
                                 const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &coeffId,
                                 UpdateType update) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto localMatrices = edge.getData(operatorId);
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto coeff = edge.getData(coeffId)->getPointer( Level );

  ValueType tmp;

  std::array<DirVertex,3> triangleGraySW = { VERTEX_C, VERTEX_W,  VERTEX_S  };
  std::array<DirVertex,3> triangleBlueS  = { VERTEX_C, VERTEX_S,  VERTEX_SE };
  std::array<DirVertex,3> triangleGraySE = { VERTEX_C, VERTEX_SE, VERTEX_E  };
  std::array<DirVertex,3> triangleGrayNW = { VERTEX_C, VERTEX_W,  VERTEX_NW };
  std::array<DirVertex,3> triangleBlueN  = { VERTEX_C, VERTEX_NW, VERTEX_N  };
  std::array<DirVertex,3> triangleGrayNE = { VERTEX_C, VERTEX_N,  VERTEX_E  };

  Face* face = storage->getFace(edge.neighborFaces()[0]);
  uint_t s_south = face->vertex_index(edge.neighborVertices()[0]);
  uint_t e_south = face->vertex_index(edge.neighborVertices()[1]);
  uint_t o_south = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));

  uint_t s_north, e_north, o_north;

  if (edge.getNumNeighborFaces() == 2) {
    face = storage->getFace(edge.neighborFaces()[1]);
    s_north = face->vertex_index(edge.neighborVertices()[0]);
    e_north = face->vertex_index(edge.neighborVertices()[1]);
    o_north = face->vertex_index(face->get_vertex_opposite_to_edge(edge.getID()));
  }

  for (size_t i = 1; i < rowsize - 1; ++i) {

    if (update == Replace) {
      tmp = ValueType(0);
    }
    else {
      tmp = dst[index<Level>(i, VERTEX_C)];
    }

    tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 0), src, coeff, triangleGraySW, {e_south, s_south, o_south});
    tmp += assembleLocal<ValueType, Level>(i, localMatrices->getBlueMatrix(Level, 0), src, coeff, triangleBlueS, {o_south, e_south, s_south});
    tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 0), src, coeff, triangleGraySE, {s_south, o_south, e_south});

    if (edge.getNumNeighborFaces() == 2) {

      tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 1), src, coeff, triangleGrayNW, {e_north, s_north, o_north});
      tmp += assembleLocal<ValueType, Level>(i, localMatrices->getBlueMatrix(Level, 1), src, coeff, triangleBlueN, {o_north, e_north, s_north});
      tmp += assembleLocal<ValueType, Level>(i, localMatrices->getGrayMatrix(Level, 1), src, coeff, triangleGrayNE, {s_north, o_north, e_north});
    }

    dst[index<Level>(i, VERTEX_C)] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, applyCoefficientTmpl, applyCoefficient )

template< typename ValueType, uint_t Level >
inline void smoothGSTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                      const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId,
                      const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &rhsId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto rhs = edge.getData(rhsId)->getPointer( Level );

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

SPECIALIZE_WITH_VALUETYPE( void, smoothGSTmpl, smooth_gs )

template< typename ValueType, uint_t Level >
inline void smoothSORTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                          const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId,
                          const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &rhsId,
                          ValueType relax) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto rhs = edge.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (size_t i = 1; i < rowsize - 1; ++i) {

    tmp = rhs[index<Level>(i, VERTEX_C)];

    for (auto& neighbor : neighbors_on_edge) {
      tmp -= opr_data[neighbor] * dst[index<Level>(i, neighbor)];
    }

    for (auto& neighbor : neighbors_south) {
      tmp -= opr_data[neighbor] * dst[index<Level>(i, neighbor)];
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto& neighbor : neighbors_north) {
        tmp -= opr_data[neighbor] * dst[index<Level>(i, neighbor)];
      }
    }
    
    dst[index<Level>(i, VERTEX_C)] = (1.0-relax) * dst[index<Level>(i, VERTEX_C)] + relax * tmp/opr_data[VERTEX_C];
  }
}

SPECIALIZE_WITH_VALUETYPE( void, smoothSORTmpl, smooth_sor )

template< typename ValueType, uint_t Level >
inline void smoothJacTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                          const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId,
                          const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &rhsId,
                          const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &tmpId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto dst = edge.getData(dstId)->getPointer( Level );
  auto rhs = edge.getData(rhsId)->getPointer( Level );
  auto tmp = edge.getData(tmpId)->getPointer( Level );

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

SPECIALIZE_WITH_VALUETYPE( void, smoothJacTmpl, smooth_jac )

template< typename ValueType, uint_t SourceLevel >
inline void prolongateTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &memoryId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize_c = levelinfo::num_microvertices_per_edge(SourceLevel);

  auto edge_data_f = edge.getData(memoryId)->getPointer( SourceLevel + 1 );
  auto edge_data_c = edge.getData(memoryId)->getPointer( SourceLevel );

  size_t i_c;
  for (i_c = 1; i_c < rowsize_c - 1; ++i_c) {

    edge_data_f[index<SourceLevel+1>(2*i_c, VERTEX_C)] = edge_data_c[index<SourceLevel>(i_c, VERTEX_C)];
    edge_data_f[index<SourceLevel+1>(2*i_c-1, VERTEX_C)] = 0.5*(edge_data_c[index<SourceLevel>(i_c-1, VERTEX_C)]
                                                              + edge_data_c[index<SourceLevel>(i_c, VERTEX_C)]);
  }

  edge_data_f[index<SourceLevel+1>(2*i_c-1, VERTEX_C)] = 0.5*(edge_data_c[index<SourceLevel>(i_c-1, VERTEX_C)]
                                                            + edge_data_c[index<SourceLevel>(i_c, VERTEX_C)]);
}

SPECIALIZE_WITH_VALUETYPE( void, prolongateTmpl, prolongate )

template< typename ValueType, uint_t Level >
inline void prolongateQuadraticTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &memoryId) {

  //TODO: rewrite using index function possible? maybe more generalized notion of Operator between different levels
  // is required.

  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(Level);
  size_t i_fine = 1;
  ValueType invtemp = 1/8.;
  const ValueType s1[3] = {3*invtemp, 6*invtemp, -invtemp};
  const ValueType s2[3] = {-invtemp, 6*invtemp, 3*invtemp};

  auto edge_data_f = edge.getData(memoryId)->getPointer( Level + 1 );
  auto edge_data_c = edge.getData(memoryId)->getPointer( Level );
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

SPECIALIZE_WITH_VALUETYPE( void, prolongateQuadraticTmpl, prolongateQuadratic )

template< typename ValueType, uint_t Level >
inline void restrictTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &memoryId) {
  using namespace EdgeCoordsVertex;

  size_t rowsize_c = levelinfo::num_microvertices_per_edge(Level - 1);

  auto edge_data_f = edge.getData(memoryId)->getPointer( Level );
  auto edge_data_c = edge.getData(memoryId)->getPointer( Level - 1 );

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

SPECIALIZE_WITH_VALUETYPE( void, restrictTmpl, restrict )

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Edge &edge, const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &dstId, uint_t& num) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for (size_t i = 1; i < rowsize - 1; ++i) {
    edge.getData(dstId)->getPointer( Level )[index<Level>(i, VERTEX_C)] = walberla::real_c(num++);
  }
}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )

template< typename ValueType, uint_t Level >
inline void integrateDGTmpl(Edge &edge,
                            const std::shared_ptr< PrimitiveStorage >& storage,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId) {

  using namespace EdgeCoordsVertex;
  typedef stencilDirection sD;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto rhs = edge.getData(rhsId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );

  ValueType tmp;

  Face* face = storage->getFace(edge.neighborFaces()[0]);
  real_t weightedFaceArea0, weightedFaceArea1;

  weightedFaceArea0 = std::pow(4.0, -walberla::real_c(Level)) * face->area / 3.0;

  uint_t s_north, e_north, o_north;

  if (edge.getNumNeighborFaces() == 2) {
    face = storage->getFace(edge.neighborFaces()[1]);
    weightedFaceArea1 = std::pow(4.0, -walberla::real_c(Level)) * face->area / 3.0;
  }

  for (size_t i = 1; i < rowsize - 1; ++i) {

    tmp =  weightedFaceArea0 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_SW)];
    tmp += weightedFaceArea0 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_BLUE_SE)];
    tmp += weightedFaceArea0 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_SE)];

    if (edge.getNumNeighborFaces() == 2) {

      tmp += weightedFaceArea1 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_NW)];
      tmp += weightedFaceArea1 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_BLUE_NW)];
      tmp += weightedFaceArea1 * rhs[DGEdge::indexDGFaceFromVertex<Level>(i, sD::CELL_GRAY_NE)];
    }

    dst[index<Level>(i, VERTEX_C)] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, integrateDGTmpl, integrateDG )

#ifdef HHG_BUILD_WITH_PETSC
template<uint_t Level>
inline void saveOperatorTmpl(Edge &edge, const PrimitiveDataID<EdgeP1StencilMemory, Edge> &operatorId,
                         const PrimitiveDataID<EdgeP1FunctionMemory< PetscInt >, Edge> &srcId,
                         const PrimitiveDataID<EdgeP1FunctionMemory< PetscInt >, Edge> &dstId, Mat& mat) {
  using namespace EdgeCoordsVertex;

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  auto &opr_data = edge.getData(operatorId)->data[Level];
  auto src = edge.getData(srcId)->getPointer( Level );
  auto dst = edge.getData(dstId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 1; ++i) {
    PetscInt dstint = dst[index<Level>(i, VERTEX_C)];
    PetscInt srcint = src[index<Level>(i, VERTEX_C)];
    //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, VERTEX_C)], opr_data[VERTEX_C]);
    MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[VERTEX_C] ,INSERT_VALUES);         //TODO: Make this more efficient by grouping all of them in an array

    for (auto& neighbor : neighbors_on_edge) {
      srcint = src[index<Level>(i, neighbor)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[neighbor] ,INSERT_VALUES);
    }

    for (auto& neighbor : neighbors_south) {
      srcint = src[index<Level>(i, neighbor)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
      MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[neighbor] ,INSERT_VALUES);
    }

    if (edge.getNumNeighborFaces() == 2) {
      for (auto& neighbor : neighbors_north) {
        srcint = src[index<Level>(i, neighbor)];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, VERTEX_C)], src[index<Level>(i, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstint,1,&srcint,&opr_data[neighbor] ,INSERT_VALUES);
      }
    }
  }
}

SPECIALIZE(void, saveOperatorTmpl, saveOperator)

template< typename ValueType, uint_t Level >
inline void createVectorFromFunctionTmpl(Edge &edge,
                                     const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &srcId,
                                     const PrimitiveDataID<EdgeP1FunctionMemory< PetscInt >, Edge> &numeratorId,
                                     Vec& vec) {
  PetscInt rowsize = (PetscInt) levelinfo::num_microvertices_per_edge(Level);

  auto src = edge.getData(srcId)->getPointer( Level );
  auto numerator = edge.getData(numeratorId)->getPointer( Level );

  VecSetValues(vec,rowsize-2,&numerator[1],&src[1],INSERT_VALUES);
}

SPECIALIZE_WITH_VALUETYPE(void, createVectorFromFunctionTmpl, createVectorFromFunction)

template< typename ValueType, uint_t Level >
inline void createFunctionFromVectorTmpl(Edge &edge,
                                         const PrimitiveDataID<EdgeP1FunctionMemory< ValueType >, Edge> &srcId,
                                         const PrimitiveDataID<EdgeP1FunctionMemory< PetscInt >, Edge> &numeratorId,
                                         Vec& vec) {
  PetscInt rowsize = (PetscInt) levelinfo::num_microvertices_per_edge(Level);

  auto numerator = edge.getData(numeratorId)->getPointer( Level );

  VecGetValues(vec,rowsize-2,&numerator[1],&edge.getData(srcId)->getPointer( Level )[1]);
}

SPECIALIZE_WITH_VALUETYPE(void, createFunctionFromVectorTmpl, createFunctionFromVector)

template< uint_t Level >
inline void applyDirichletBCTmpl(Edge &edge,std::vector<PetscInt> &mat,
                                 const PrimitiveDataID<EdgeP1FunctionMemory< PetscInt >, Edge> &numeratorId){

  size_t rowsize = levelinfo::num_microvertices_per_edge(Level);

  for(uint_t i = 1;i<rowsize-1; i++)
  {
    mat.push_back(edge.getData(numeratorId)->getPointer( Level )[i]);
  }

}
SPECIALIZE(void, applyDirichletBCTmpl, applyDirichletBC)
#endif

}
}
