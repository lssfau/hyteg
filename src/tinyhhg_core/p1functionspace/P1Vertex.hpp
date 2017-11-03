
#pragma once

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

namespace hhg {

namespace P1Vertex {

template< typename ValueType >
inline void interpolate(Vertex &vertex,
                        const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &vertexMemoryId,
                        const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &srcIds,
                        std::function<ValueType(const hhg::Point3D &, const std::vector<ValueType>&)> &expr,
                        uint_t level) {
  VertexP1FunctionMemory< ValueType > *vertexMemory = vertex.getData(vertexMemoryId);
  std::vector<ValueType> srcVector(srcIds.size());

  for (uint_t k = 0; k < srcIds.size(); ++k) {
    srcVector[k] = vertex.getData(srcIds[k])->getPointer(level)[0];
  }

  vertexMemory->getPointer( level )[0] = expr(vertex.getCoordinates(), srcVector);
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
inline void smooth_sor(Vertex &vertex, const PrimitiveDataID<VertexP1StencilMemory, Vertex> &operatorId,
                      const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &dstId,
                      const PrimitiveDataID<VertexP1FunctionMemory< ValueType >, Vertex> &rhsId, size_t level,
                      ValueType relax) {
  auto &opr_data = vertex.getData(operatorId)->data[level];
  auto dst = vertex.getData(dstId)->getPointer( level );
  auto rhs = vertex.getData(rhsId)->getPointer( level );

  ValueType tmp;
  tmp = rhs[0];

  for (size_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
    tmp -= opr_data[i + 1]*dst[i + 1];
  }

  dst[0] = (1.0-relax) * dst[0] + relax * tmp/opr_data[0];
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
inline void enumerate(size_t level, Vertex &vertex, const PrimitiveDataID <VertexP1FunctionMemory<ValueType>, Vertex> &dstId, uint_t &num) {
  auto dst = vertex.getData(dstId)->getPointer( level );
  dst[0] = static_cast< ValueType >( num++ );
}

template< typename ValueType >
inline void integrateDG(Vertex &vertex,
                            const std::shared_ptr< PrimitiveStorage >& storage,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &rhsP1Id,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                            uint_t level) {

  auto rhs = vertex.getData(rhsId)->getPointer( level );
  auto rhsP1 = vertex.getData(rhsP1Id)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );

  ValueType tmp = 0.0;

  for(auto faceIt : vertex.neighborFaces()) {
    Face *face = storage->getFace(faceIt.getID());

    real_t weightedFaceArea = std::pow(4.0, -walberla::real_c(level))*face->area / 3.0;

    uint_t localFaceId = vertex.face_index(face->getID());

    uint_t faceMemoryIndex = 2 * localFaceId;

    std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());
    uint_t edge_idx[2] = { vertex.edge_index(adj_edges[0]) + 1, vertex.edge_index(adj_edges[1]) + 1 };

    tmp += weightedFaceArea * rhs[faceMemoryIndex] * (0.5 * 0.5 * (rhsP1[0] + rhsP1[edge_idx[0]]) + 0.5 * 0.5 * (rhsP1[0] + rhsP1[edge_idx[1]]));
  }

  dst[0] = tmp;
}

template< typename ValueType >
inline void projectNormal(Vertex &vertex,
                          const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &velocityIds,
                          const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &normalIds,
                          uint_t level) {

  Point2D normal({ vertex.getData(normalIds[0])->getPointer(level)[0], vertex.getData(normalIds[1])->getPointer(level)[0] });
  Point2D velocity({ vertex.getData(velocityIds[0])->getPointer(level)[0], vertex.getData(velocityIds[1])->getPointer(level)[0] });

  Matrix2r projector;
  projector(0,0) = 1.0 - normal[0]*normal[0];
  projector(0,1) = -normal[0]*normal[1];
  projector(1,0) = -normal[1]*normal[0];
  projector(1,1) = 1.0 - normal[1]*normal[1];

  Point2D projected = projector.multiply(velocity);

  vertex.getData(velocityIds[0])->getPointer(level)[0] = projected(0);
  vertex.getData(velocityIds[1])->getPointer(level)[0] = projected(1);
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
