
#pragma once

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

namespace hhg {
namespace edgedof {
namespace macroedge {

using walberla::uint_t;
using walberla::real_c;

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Edge & edge,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & edgeMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>> &srcIds,
                            std::function< ValueType( const hhg::Point3D &, const std::vector<ValueType>& ) > & expr)
{
  auto edgeData = edge.getData( edgeMemoryId )->getPointer( Level );

  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(edge.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVector(srcIds.size());

  const Point3D leftCoords  = edge.getCoordinates()[0];
  const Point3D rightCoords = edge.getCoordinates()[1];

  const Point3D microEdgeOffset = ( rightCoords - leftCoords ) / real_c( 2 * levelinfo::num_microedges_per_edge( Level ) );

  for ( const auto & it : edgedof::macroedge::Iterator( Level ) )
  {
    const Point3D currentCoordinates = leftCoords + microEdgeOffset + 2 * it.col() * microEdgeOffset;

    for (uint_t k = 0; k < srcPtr.size(); ++k) {
      srcVector[k] = srcPtr[k][edgedof::macroedge::horizontalIndex( Level, it.col())];
    }

    edgeData[edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C )] = expr( currentCoordinates, srcVector );
  }
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate )


template< typename ValueType, uint_t Level >
inline void addTmpl( Edge & edge, const std::vector< ValueType > & scalars,
                     const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > & srcIds,
                     const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );
  WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" );

  auto dstData = edge.getData( dstId )->getPointer( Level );

  for ( const auto & it : edgedof::macroedge::Iterator( Level ) )
  {
    ValueType tmp = static_cast< ValueType >( 0.0 );

    const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C );

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = edge.getData( srcIds[i] )->getPointer( Level );

      tmp += scalar * srcData[ idx ];
    }

    dstData[ idx ] += tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add )


template< typename ValueType, uint_t Level >
inline void assignTmpl( Edge & edge, const std::vector< ValueType > & scalars,
                        const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > & srcIds,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Edge > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );
  WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" );

  auto dstData = edge.getData( dstId )->getPointer( Level );

  for ( const auto & it : edgedof::macroedge::Iterator( Level ) )
  {
    ValueType tmp = static_cast< ValueType >( 0.0 );

    const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C );

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = edge.getData( srcIds[i] )->getPointer( Level );

      tmp += scalar * srcData[ idx ];
    }

    dstData[ idx ] = tmp;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, assignTmpl, assign )


template< typename ValueType, uint_t Level >
inline real_t dotTmpl( Edge & edge,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& lhsId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& rhsId )
{
  auto lhsData = edge.getData( lhsId )->getPointer( Level );
  auto rhsData = edge.getData( rhsId )->getPointer( Level );

  real_t scalarProduct = real_c( 0 );

  for ( const auto & it : edgedof::macroedge::Iterator( Level ) )
  {
    const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C );
    scalarProduct += lhsData[ idx ] * rhsData[ idx ];
  }

  return scalarProduct;
}

SPECIALIZE_WITH_VALUETYPE( real_t, dotTmpl, dot )


template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Edge &edge,
                          const PrimitiveDataID < FunctionMemory< ValueType >, Edge> &dstId,
                          uint_t& num)
{
  ValueType *dst = edge.getData(dstId)->getPointer(Level);

  for(uint_t i = 0 ; i < levelinfo::num_microedges_per_edge( Level ) ; ++i){
    dst[hhg::edgedof::macroedge::horizontalIndex( Level, i )] = num;
    ++num;
  }
}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )

template<uint_t Level>
inline void applyTmpl(Edge &edge,
                       const PrimitiveDataID<StencilMemory < real_t >, Edge> &operatorId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &srcId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Edge> &dstId,
                       UpdateType update)
{
  using namespace hhg::edgedof::macroedge;
  size_t rowsize = levelinfo::num_microedges_per_edge(Level);

  real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  real_t * src      = edge.getData(srcId)->getPointer( Level );
  real_t * dst      = edge.getData(dstId)->getPointer( Level );

  real_t tmp;

  for(uint_t i = 0; i < rowsize; ++i){
    tmp = 0.0;
    for(uint_t k = 0; k < neighborsOnEdgeFromHorizontalEdge.size(); ++k){
      tmp += opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnEdgeFromHorizontalEdge[k])] *
             src[indexFromHorizontalEdge( Level, i, neighborsOnEdgeFromHorizontalEdge[k] )];
    }
    for(uint_t k = 0; k < neighborsOnSouthFaceFromHorizontalEdge.size(); ++k){
      tmp += opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnSouthFaceFromHorizontalEdge[k])] *
             src[indexFromHorizontalEdge( Level, i, neighborsOnSouthFaceFromHorizontalEdge[k] )];
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < neighborsOnNorthFaceFromHorizontalEdge.size(); ++k){
        tmp += opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnNorthFaceFromHorizontalEdge[k])] *
               src[indexFromHorizontalEdge( Level, i, neighborsOnNorthFaceFromHorizontalEdge[k] )];
      }
    }

    if (update==Replace) {
      dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] = tmp;
    } else if (update==Add) {
      dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )] += tmp;
    }
  }
}

SPECIALIZE(void, applyTmpl, apply)

template< typename ValueType, size_t Level >
inline void printFunctionMemory(Edge& edge, const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &dstId){
  ValueType* edgeMemory = edge.getData(dstId)->getPointer( Level );
  using namespace std;
  cout << setfill('=') << setw(100) << "" << endl;
  cout << edge << std::left << setprecision(1) << fixed << setfill(' ') << endl;
  uint_t rowsize = levelinfo::num_microvertices_per_edge( Level );
  cout << "Horizontal Edge" << endl;
  if(edge.getNumNeighborFaces() == 2) {
    for (uint_t i = 1; i < rowsize-1; ++i) {
      cout << setw(5) << edgeMemory[hhg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_NW )] << "|";
    }
    cout << endl;
  }
  for(uint_t i = 1; i < rowsize; ++i){
    cout << setw(5) << edgeMemory[hhg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_W )] << "|";
  }
  cout << endl << "     |";
  for(uint_t i = 1; i < rowsize-1; ++i){
    cout << setw(5) << edgeMemory[hhg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_HO_SE )] << "|";
  }
  cout << endl << "Diagonal Edge" << endl;
  if(edge.getNumNeighborFaces() == 2) {
    for (uint_t i = 1; i < rowsize; ++i) {
      cout << setw(5) << edgeMemory[hhg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_DI_NW )] << "|";
    }
    cout << endl;
  }
  for(uint_t i = 0; i < rowsize-1; ++i){
    cout << setw(5) << edgeMemory[hhg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_DI_SE )] << "|";
  }
  cout << endl << "Vertical Edge" << endl;
  if(edge.getNumNeighborFaces() == 2) {
    for (uint_t i = 0; i < rowsize -1; ++i) {
      cout << setw(5) << edgeMemory[hhg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_VE_N )] << "|";
    }
    cout << endl;
  }
  for(uint_t i = 1; i < rowsize; ++i){
    cout << setw(5) << edgeMemory[hhg::edgedof::macroedge::indexFromVertex( Level, i, stencilDirection::EDGE_VE_S )] << "|";
  }
  cout << endl << setfill('=') << setw(100) << "" << endl << setfill(' ');

}

#ifdef HHG_BUILD_WITH_PETSC
template< typename ValueType, uint_t Level >
inline void createVectorFromFunctionTmpl(Edge &edge,
                                         const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId,
                                         Vec& vec) {
  auto src = edge.getData(srcId)->getPointer( Level );
  auto numerator = edge.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : edgedof::macroedge::Iterator( Level ) )
  {
    const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C );
    VecSetValues(vec,1,&numerator[idx],&src[idx],INSERT_VALUES);
  }
}

SPECIALIZE_WITH_VALUETYPE(void, createVectorFromFunctionTmpl, createVectorFromFunction)

template< typename ValueType, uint_t Level >
inline void createFunctionFromVectorTmpl(Edge &edge,
                                         const PrimitiveDataID<FunctionMemory< ValueType >, Edge> &srcId,
                                         const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId,
                                         Vec& vec) {
  auto src = edge.getData(srcId)->getPointer( Level );
  auto numerator = edge.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : edgedof::macroedge::Iterator( Level ) )
  {
    const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C );
    VecGetValues(vec,1,&numerator[idx],&src[idx]);
  }

}

SPECIALIZE_WITH_VALUETYPE(void, createFunctionFromVectorTmpl, createFunctionFromVector)

template< uint_t Level >
inline void applyDirichletBCTmpl(Edge &edge,std::vector<PetscInt> &mat,
                                 const PrimitiveDataID<FunctionMemory< PetscInt >, Edge> &numeratorId){

  auto numerator = edge.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : edgedof::macroedge::Iterator( Level ) )
  {
    const uint_t idx = edgedof::macroedge::indexFromHorizontalEdge( Level, it.col(), stencilDirection::EDGE_HO_C );
    mat.push_back(numerator[idx]);
  }

}
SPECIALIZE(void, applyDirichletBCTmpl, applyDirichletBC)
#endif


} ///namespace macroedge
} ///namespace edgedof
} ///namespace hhg
