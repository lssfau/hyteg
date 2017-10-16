
#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

namespace hhg {

namespace P1Vertex {

template< typename ValueType >
inline void interpolate(Vertex &vertex,
                        const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &vertexMemoryId,
                        std::function<ValueType(const hhg::Point3D &)> &expr,
                        size_t level) {
  VertexP1FunctionMemory< ValueType > *vertexMemory = vertex.getData(vertexMemoryId);
  vertexMemory->getPointer( level )[0] = expr(vertex.getCoordinates());
}

template< typename ValueType >
inline void assign(Vertex &vertex,
                   const std::vector<ValueType> &scalars,
                   const std::vector<PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex>> &srcIds,
                   const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId,
                   size_t level) {
  ValueType tmp = scalars[0]*vertex.getData(srcIds[0])->getPointer( level )[0];

  for (size_t i = 1; i < srcIds.size(); ++i) {
    tmp += scalars[i]*vertex.getData(srcIds[i])->getPointer( level )[0];
  }

  vertex.getData(dstId)->getPointer( level )[0] = tmp;
}

template< typename ValueType >
inline void add(Vertex &vertex,
                const std::vector<ValueType> &scalars,
                const std::vector<PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex>> &srcIds,
                const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId,
                size_t level) {
  ValueType tmp = 0.0;

  for (size_t i = 0; i < srcIds.size(); ++i) {
    tmp += scalars[i]*vertex.getData(srcIds[i])->getPointer( level )[0];
  }

  vertex.getData(dstId)->getPointer( level )[0] += tmp;
}

template< typename ValueType >
inline real_t dot(Vertex &vertex,
                  const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &lhsMemoryId,
                  const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &rhsMemoryId,
                  size_t level) {
  return vertex.getData(lhsMemoryId)->getPointer( level )[0]*vertex.getData(rhsMemoryId)->getPointer( level )[0];
}

template< typename ValueType >
inline void apply(Vertex &vertex,
                  const PrimitiveDataID<VertexP1StencilMemory, Vertex> &operatorId,
                  const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &srcId,
                  const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId,
                  size_t level,
                  UpdateType update) {
  auto &opr_data = vertex.getData(operatorId)->data[level];
  auto src = vertex.getData(srcId)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );

  if (update==Replace) {
    dst[0] = opr_data[0]*src[0];
  } else if (update==Add) {
    dst[0] += opr_data[0]*src[0];
  }

  for (size_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
    dst[0] += opr_data[i + 1]*src[i + 1];
  }
}

/// Apply function in the case of a coefficient
template< typename ValueType >
inline void applyCoefficient(Vertex &vertex,
                             const std::shared_ptr< PrimitiveStorage >& storage,
                             const PrimitiveDataID<VertexP1LocalMatrixMemory, Vertex> &operatorId,
                             const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &srcId,
                             const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId,
                             const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &coeffId,
                             uint_t level,
                             UpdateType update) {
  auto localMatrices = vertex.getData(operatorId);
  auto src = vertex.getData(srcId)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );
  auto coeff = vertex.getData(coeffId)->getPointer( level );

  if (update == Replace) {
    dst[0] = real_t(0);
  }

  uint_t neighborId = 0;
  for (auto& faceId : vertex.neighborFaces()) {
    real_t meanCoefficient = coeff[0];
    Matrix3r& local_stiffness = localMatrices->getGrayMatrix(level, neighborId);

    Face* face = storage->getFace(faceId);

    uint_t v_i = face->vertex_index(vertex.getID());

    std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());

    for (auto &edgeId : adj_edges) {
      uint_t edge_idx = vertex.edge_index(edgeId) + 1;
      meanCoefficient += coeff[edge_idx];
    }

    meanCoefficient *= 1.0/3.0;

    // iterate over adjacent edges
    for (auto &edgeId : adj_edges) {
      uint_t edge_idx = vertex.edge_index(edgeId) + 1;
      Edge *edge = storage->getEdge(edgeId);
      PrimitiveID vertex_j = edge->get_opposite_vertex(vertex.getID());

      uint_t v_j = face->vertex_index(vertex_j);

      dst[0] += meanCoefficient * local_stiffness(v_i, v_j) * src[edge_idx];
    }

    // add contribution of center vertex
    dst[0] += meanCoefficient * local_stiffness(v_i, v_i) * src[0];
    ++neighborId;
  }
}

template< typename ValueType >
inline void smooth_gs(Vertex &vertex, const PrimitiveDataID<VertexP1StencilMemory, Vertex> &operatorId,
                      const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId,
                      const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &rhsId, size_t level) {
  auto &opr_data = vertex.getData(operatorId)->data[level];
  auto dst = vertex.getData(dstId)->getPointer( level );
  auto rhs = vertex.getData(rhsId)->getPointer( level );

  dst[0] = rhs[0];

  for (size_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
    dst[0] -= opr_data[i + 1]*dst[i + 1];
  }

  dst[0] /= opr_data[0];
}

template< typename ValueType >
inline void smooth_jac(Vertex &vertex, const PrimitiveDataID<VertexP1StencilMemory, Vertex> &operatorId,
                      const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId,
                      const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &rhsId,
                      const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &tmpId, size_t level) {
  auto &opr_data = vertex.getData(operatorId)->data[level];
  auto dst = vertex.getData(dstId)->getPointer( level );
  auto rhs = vertex.getData(rhsId)->getPointer( level );
  auto tmp = vertex.getData(tmpId)->getPointer( level );

  dst[0] = rhs[0];

  for (size_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
    dst[0] -= opr_data[i + 1]*tmp[i + 1];
  }

  dst[0] /= opr_data[0];
}

template< typename ValueType >
inline void prolongate(Vertex &vertex, const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &memoryId, size_t sourceLevel) {
  vertex.getData(memoryId)->getPointer(sourceLevel + 1)[0] =
      vertex.getData(memoryId)->getPointer(sourceLevel)[0];
}

template< typename ValueType >
inline void prolongateQuadratic(Vertex &vertex,
                                const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &memoryId,
                                size_t level) {
  prolongate(vertex, memoryId, level);
}

template< typename ValueType >
inline void restrict(Vertex &vertex, const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &memoryId, size_t level) {
  auto vertex_data_f = vertex.getData(memoryId)->getPointer( level );
  auto vertex_data_c = vertex.getData(memoryId)->getPointer( level - 1 );

  vertex_data_c[0] = vertex_data_f[0];

  for (uint_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
    vertex_data_c[0] += 0.5*vertex_data_f[i+1];
    i += 1;
  }
}

template< typename ValueType >
inline void enumerate(Vertex &vertex, const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId, size_t level, uint_t& num) {
  auto dst = vertex.getData(dstId)->getPointer( level );
  dst[0] = static_cast< ValueType >( num++ );
}

#ifdef HHG_BUILD_WITH_PETSC
inline void saveOperator(Vertex &vertex,
                         const PrimitiveDataID<VertexP1StencilMemory, Vertex> &operatorId,
                         const PrimitiveDataID<VertexP1FunctionMemory< PetscInt >, Vertex> &srcId,
                         const PrimitiveDataID<VertexP1FunctionMemory< PetscInt >, Vertex> &dstId,
                         Mat& mat,
                         uint_t level) {
  auto &opr_data = vertex.getData(operatorId)->data[level];
  auto src = vertex.getData(srcId)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );

  MatSetValues(mat, 1, dst, (PetscInt) (vertex.getNumNeighborEdges() + 1), src, opr_data.get(), INSERT_VALUES);
}

template< typename ValueType >
inline void createVectorFromFunction(Vertex &vertex,
                         const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &srcId,
                         const PrimitiveDataID<VertexP1FunctionMemory< PetscInt >, Vertex> &numeratorId,
                         Vec& vec,
                         uint_t level) {

  auto src = vertex.getData(srcId)->getPointer( level );
  PetscInt numerator = vertex.getData(numeratorId)->getPointer( level )[0];

  VecSetValues(vec, 1, &numerator, src, INSERT_VALUES);

}

template< typename ValueType >
inline void createFunctionFromVector(Vertex &vertex,
                                     const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &srcId,
                                     const PrimitiveDataID<VertexP1FunctionMemory< PetscInt >, Vertex> &numeratorId,
                                     Vec& vec,
                                     uint_t level) {


  PetscInt numerator = vertex.getData(numeratorId)->getPointer( level )[0];

  VecGetValues(vec, 1, &numerator, vertex.getData(srcId)->getPointer( level ));

}

inline void applyDirichletBC(Vertex &vertex,std::vector<PetscInt> &mat, uint_t level,
                             const PrimitiveDataID<VertexP1FunctionMemory< PetscInt >, Vertex> &numeratorId){

  mat.push_back(vertex.getData(numeratorId)->getPointer( level )[0]);

}

#endif


}
}
