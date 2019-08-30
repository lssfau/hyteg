
#pragma once

#include "hyteg/Levelinfo.hpp"

#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"

#ifdef DEBUG_ELEMENTWISE
#include "hyteg/format.hpp"
#endif

namespace hyteg {
namespace vertexdof {
namespace macrovertex {

template< typename ValueType >
inline void interpolate( const uint_t & level, const Vertex & vertex,
                         const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &vertexMemoryId,
                         const ValueType & scalar )
{
  auto vertexMemory = vertex.getData(vertexMemoryId)->getPointer(level);
  vertexMemory[0] = scalar;
}


template< typename ValueType >
inline void interpolate(Vertex &vertex,
                        const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &vertexMemoryId,
                        const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &srcIds,
                        const std::function<ValueType(const hyteg::Point3D &, const std::vector<ValueType>&)> &expr,
                        uint_t level) {
  FunctionMemory< ValueType > *vertexMemory = vertex.getData(vertexMemoryId);
  std::vector<ValueType> srcVector(srcIds.size());

  for (uint_t k = 0; k < srcIds.size(); ++k) {
    srcVector[k] = vertex.getData(srcIds[k])->getPointer(level)[0];
  }

  Point3D xBlend;
  vertex.getGeometryMap()->evalF(vertex.getCoordinates(), xBlend);
  vertexMemory->getPointer( level )[0] = expr(xBlend, srcVector);
}

template< typename ValueType >
inline void swap( const uint_t & level, Vertex & vertex,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex > & dstID )
{
  auto srcData = vertex.getData( srcID );
  auto dstData = vertex.getData( dstID );
  srcData->swap( *dstData, level );
}

template< typename ValueType >
inline void assign(Vertex &vertex,
                   const std::vector<ValueType> &scalars,
                   const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &srcIds,
                   const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                   size_t level) {
  ValueType tmp = scalars[0]*vertex.getData(srcIds[0])->getPointer( level )[0];

  for (size_t i = 1; i < srcIds.size(); ++i) {
    tmp += scalars[i]*vertex.getData(srcIds[i])->getPointer( level )[0];
  }

  vertex.getData(dstId)->getPointer( level )[0] = tmp;
}

template< typename ValueType >
inline void add(const Vertex & vertex,
                const ValueType & scalar,
                const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> & dstId,
                const uint_t & level)
{
  vertex.getData(dstId)->getPointer( level )[0] += scalar;
}


template< typename ValueType >
inline void add(Vertex &vertex,
                const std::vector<ValueType> &scalars,
                const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &srcIds,
                const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                size_t level) {
  auto tmp = ValueType( 0 );

  for (size_t i = 0; i < srcIds.size(); ++i) {
    tmp += scalars[i]*vertex.getData(srcIds[i])->getPointer( level )[0];
  }

  vertex.getData(dstId)->getPointer( level )[0] += tmp;
}

template< typename ValueType >
inline void multElementwise(Vertex &vertex,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> &srcIds,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                            size_t level) {
  ValueType tmp = vertex.getData(srcIds[0])->getPointer( level )[0];

  for (size_t i = 1; i < srcIds.size(); ++i) {
    tmp *= vertex.getData(srcIds[i])->getPointer( level )[0];
  }

  vertex.getData(dstId)->getPointer( level )[0] = tmp;
}

template< typename ValueType >
inline ValueType dot(Vertex &vertex,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &lhsMemoryId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &rhsMemoryId,
                  size_t level) {
  return vertex.getData(lhsMemoryId)->getPointer( level )[0]*vertex.getData(rhsMemoryId)->getPointer( level )[0];
}

template < typename ValueType >
inline ValueType sum( const uint_t&                                                 level,
                      const Vertex&                                                 vertex,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& dataID,
                      const bool&                                                   absolute )
{
   if ( absolute )
   {
      return std::abs( vertex.getData( dataID )->getPointer( level )[0] );
   }
   return vertex.getData( dataID )->getPointer( level )[0];
}

template< typename ValueType >
inline void apply(Vertex &vertex,
                  const PrimitiveDataID<StencilMemory< ValueType >, Vertex> &operatorId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                  size_t level,
                  UpdateType update) {
  auto opr_data = vertex.getData(operatorId)->getPointer( level );
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

template< typename ValueType >
inline void applyPointwise(const uint_t & level,
                  const Vertex &vertex,
                  const PrimitiveDataID<StencilMemory< ValueType >, Vertex> &operatorId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                  UpdateType update) {
  auto opr_data = vertex.getData(operatorId)->getPointer( level );
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

template< typename ValueType >
inline void smooth_gs(Vertex &vertex, const PrimitiveDataID<StencilMemory< ValueType >, Vertex> &operatorId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &rhsId, size_t level) {
  auto opr_data = vertex.getData( operatorId )->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );
  auto rhs = vertex.getData(rhsId)->getPointer( level );

  dst[0] = rhs[0];

  for (size_t i = 0; i < vertex.getNumNeighborEdges(); ++i) {
    dst[0] -= opr_data[i + 1]*dst[i + 1];
  }

  dst[0] /= opr_data[0];
}

template< typename ValueType >
inline void smooth_sor(Vertex &vertex, const PrimitiveDataID<StencilMemory< ValueType >, Vertex> &operatorId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &rhsId, size_t level,
                      ValueType relax) {
  auto opr_data = vertex.getData(operatorId)->getPointer( level );
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
inline void smooth_jac(Vertex &vertex, const PrimitiveDataID<StencilMemory< ValueType >, Vertex> &operatorId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &dstId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &rhsId,
                      const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &tmpId, size_t level) {
  auto opr_data = vertex.getData(operatorId)->getPointer( level );
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
inline void enumerate(size_t level, Vertex &vertex, const PrimitiveDataID <FunctionMemory<ValueType>, Vertex> &dstId, ValueType &num) {
  auto dst = vertex.getData(dstId)->getPointer( level );
  dst[0] = num++ ;
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

  auto tmp = real_t( 0 );

  for(const auto& faceIt : vertex.neighborFaces()) {
    Face *face = storage->getFace(faceIt.getID());

    real_t weightedFaceArea = std::pow(4.0, -walberla::real_c(level))*face->area / 3.0;

    uint_t localFaceId = vertex.face_index(face->getID());

    uint_t faceMemoryIndex = 2 * localFaceId;

    std::vector<PrimitiveID> adj_edges = face->adjacent_edges(vertex.getID());
    uint_t edge_idx[2] = { vertex.edge_index(adj_edges[0]) + 1, vertex.edge_index(adj_edges[1]) + 1 };

    tmp += weightedFaceArea * rhs[faceMemoryIndex] * (0.5 * 0.5 * (rhsP1[0] + rhsP1[edge_idx[0]]) + 0.5 * 0.5 * (rhsP1[0] + rhsP1[edge_idx[1]]));
  }

  dst[0] = ValueType( tmp );
}


template< typename ValueType >
inline ValueType getMaxValue( const uint_t &level, Vertex &vertex, const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId ) {
  auto src = vertex.getData( srcId )->getPointer( level );
  return src[0];
}


template< typename ValueType >
inline ValueType getMaxMagnitude( const uint_t &level, Vertex &vertex, const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId ) {
  auto src = vertex.getData( srcId )->getPointer( level );
  return std::abs( src[0] );
}

template< typename ValueType >
inline ValueType getMinValue( const uint_t &level, Vertex &vertex, const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId ) {
  auto src = vertex.getData( srcId )->getPointer( level );
  return src[0];
}


#ifdef HHG_BUILD_WITH_PETSC
inline void saveOperator(Vertex &vertex,
                         const PrimitiveDataID<StencilMemory< real_t >, Vertex> &operatorId,
                         const PrimitiveDataID<FunctionMemory< PetscInt >, Vertex> &srcId,
                         const PrimitiveDataID<FunctionMemory< PetscInt >, Vertex> &dstId,
                         Mat& mat,
                         uint_t level) {
  auto opr_data = vertex.getData(operatorId)->getPointer( level );
  auto src = vertex.getData(srcId)->getPointer( level );
  auto dst = vertex.getData(dstId)->getPointer( level );

  MatSetValues(mat, 1, dst, (PetscInt) (vertex.getNumNeighborEdges() + 1), src, opr_data, ADD_VALUES);
}

template < typename ValueType >
inline void createVectorFromFunction( const Vertex&                                                 vertex,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& srcId,
                                      const PrimitiveDataID< FunctionMemory< PetscInt >, Vertex >&  numeratorId,
                                      Vec&                                                          vec,
                                      uint_t                                                        level )
{
   auto     src       = vertex.getData( srcId )->getPointer( level );
   PetscInt numerator = vertex.getData( numeratorId )->getPointer( level )[0];

   VecSetValues( vec, 1, &numerator, src, INSERT_VALUES );
}

template< typename ValueType >
inline void createFunctionFromVector(Vertex &vertex,
                                     const PrimitiveDataID<FunctionMemory< ValueType >, Vertex> &srcId,
                                     const PrimitiveDataID<FunctionMemory< PetscInt >, Vertex> &numeratorId,
                                     Vec& vec,
                                     uint_t level) {


  PetscInt numerator = vertex.getData(numeratorId)->getPointer( level )[0];

  VecGetValues(vec, 1, &numerator, vertex.getData(srcId)->getPointer( level ));

}

inline void applyDirichletBC(Vertex &vertex,std::vector<PetscInt> &mat, uint_t level,
                             const PrimitiveDataID<FunctionMemory< PetscInt >, Vertex> &numeratorId){

  mat.push_back(vertex.getData(numeratorId)->getPointer( level )[0]);

}

#endif

}
}
}
