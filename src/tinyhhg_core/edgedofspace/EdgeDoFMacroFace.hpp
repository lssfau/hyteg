
#pragma once


#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>
#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Face & face,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Face > & faceMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                            std::function< ValueType( const hhg::Point3D &, const std::vector<ValueType>& ) > & expr)
{
  auto faceData = face.getData( faceMemoryId )->getPointer( Level );

  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVectorHorizontal(srcIds.size());
  std::vector<ValueType> srcVectorVertical(srcIds.size());
  std::vector<ValueType> srcVectorDiagonal(srcIds.size());

  const Point3D faceBottomLeftCoords  = face.coords[0];
  const Point3D faceBottomRightCoords = face.coords[1];
  const Point3D faceTopLeftCoords     = face.coords[2];

  const Point3D horizontalMicroEdgeOffset = ( ( faceBottomRightCoords - faceBottomLeftCoords ) / levelinfo::num_microedges_per_edge( Level ) ) * 0.5;
  const Point3D verticalMicroEdgeOffset   = ( ( faceTopLeftCoords     - faceBottomLeftCoords ) / levelinfo::num_microedges_per_edge( Level ) ) * 0.5;

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( ( it.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( it.row() * 2     ) * verticalMicroEdgeOffset );
    const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( ( it.col() * 2     ) * horizontalMicroEdgeOffset + ( it.row() * 2 + 1 ) * verticalMicroEdgeOffset );
    const Point3D diagonalMicroEdgePosition   = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

    // Do not update horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
        srcVectorHorizontal[k] = srcPtr[k][edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() )];
      }

      faceData[ edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = expr( horizontalMicroEdgePosition, srcVectorHorizontal );
    }

    // Do not update vertical DoFs at left border
    if ( it.col() != 0 )
    {
      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
        srcVectorVertical[k] = srcPtr[k][edgedof::macroface::verticalIndex< Level >( it.col(), it.row() )];
      }

      faceData[ edgedof::macroface::verticalIndex< Level >  ( it.col(), it.row() ) ] = expr( verticalMicroEdgePosition, srcVectorVertical );
    }

    // Do not update diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
        srcVectorDiagonal[k] = srcPtr[k][edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() )];
      }

      faceData[ edgedof::macroface::diagonalIndex< Level >  ( it.col(), it.row() ) ] = expr( diagonalMicroEdgePosition, srcVectorDiagonal );
    }
  }
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate )


template< typename ValueType, uint_t Level >
inline void addTmpl( Face & face, const std::vector< ValueType > & scalars,
                     const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcIds,
                     const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );

  auto dstData = face.getData( dstId )->getPointer( Level );

  for ( uint_t i = 0; i < scalars.size(); i++ )
  {
    const real_t scalar  = scalars[i];
    auto         srcData = face.getData( srcIds[i] )->getPointer( Level );

    for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
    {
      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
        const uint_t idx = edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() );
        dstData[ idx ] += scalar * srcData[ idx ];
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
        const uint_t idx = edgedof::macroface::verticalIndex< Level >( it.col(), it.row() );
        dstData[ idx ] += scalar * srcData[ idx ];
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
        const uint_t idx = edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() );
        dstData[ idx ] += scalar * srcData[ idx ];
      }
    }
  }

}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add )

template< typename ValueType, uint_t Level >
inline void assignTmpl( Face & face, const std::vector< ValueType > & scalars,
                        const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcIds,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );

  auto dstData = face.getData( dstId )->getPointer( Level );

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    // Do not update horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      const uint_t idx = edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() );
      dstData[ idx ] = static_cast< ValueType >( 0 );
    }

    // Do not update vertical DoFs at left border
    if ( it.col() != 0 )
    {
      const uint_t idx = edgedof::macroface::verticalIndex< Level >( it.col(), it.row() );
      dstData[ idx ] = static_cast< ValueType >( 0 );
    }

    // Do not update diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      const uint_t idx = edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() );
      dstData[ idx ] = static_cast< ValueType >( 0 );
    }
  }

  addTmpl< ValueType, Level >( face, scalars, srcIds, dstId );
}

SPECIALIZE_WITH_VALUETYPE( void, assignTmpl, assign )

template< typename ValueType, uint_t Level >
inline real_t dotTmpl( Face & face,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& lhsId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId )
{
  auto lhsData = face.getData( lhsId )->getPointer( Level );
  auto rhsData = face.getData( rhsId )->getPointer( Level );

  real_t scalarProduct = real_c( 0 );

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    // Do not read horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      const uint_t idx = edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    // Do not read vertical DoFs at left border
    if ( it.col() != 0 )
    {
      const uint_t idx = edgedof::macroface::verticalIndex< Level >( it.col(), it.row() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    // Do not read diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      const uint_t idx = edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }
  }

  return scalarProduct;
}

SPECIALIZE_WITH_VALUETYPE( real_t, dotTmpl, dot )

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Face &face,
                          const PrimitiveDataID < FunctionMemory< ValueType >, Face> &dstId,
                          uint_t& num)
{
  ValueType *dst = face.getData(dstId)->getPointer(Level);
  size_t horizontal_num = num;
  size_t diagonal_num = num +
                        hhg::edgedof::levelToFaceSizeAnyEdgeDoF< Level > -
                        hhg::levelinfo::num_microedges_per_edge( Level ) ;
  size_t vertical_num = num +
                        (hhg::edgedof::levelToFaceSizeAnyEdgeDoF< Level > -
                        hhg::levelinfo::num_microedges_per_edge( Level ))  *
                        2;
  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    /// the border edge DoFs belong to the corresponding edges
    if( it.row() != 0) {
      dst[hhg::edgedof::macroface::horizontalIndex< Level >(it.col(), it.row())] = horizontal_num;
      ++horizontal_num;
      ++num;
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      dst[hhg::edgedof::macroface::diagonalIndex< Level >(it.col(), it.row())] = diagonal_num;
      ++diagonal_num;
      ++num;
    }
    if( it.col() != 0) {
      dst[hhg::edgedof::macroface::verticalIndex< Level >(it.col(), it.row())] = vertical_num;
      ++vertical_num;
      ++num;
    }

  }

}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate )

template<uint_t Level>
inline void applyTmpl(Face &face,
                       const PrimitiveDataID<StencilMemory < real_t >, Face> &operatorId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &srcId,
                       const PrimitiveDataID<FunctionMemory< real_t >, Face> &dstId,
                       UpdateType update)
{
  real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  real_t * src      = face.getData(srcId)->getPointer( Level );
  real_t * dst      = face.getData(dstId)->getPointer( Level );

  real_t tmp;

  using namespace edgedof::macroface;

  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromHorizontalEdge.size(); ++k){
        tmp += opr_data[edgedof::stencilIndexFromHorizontalEdge(neighborsFromHorizontalEdge[k])] *
               src[indexFromHorizontalEdge< Level >(it.col(), it.row(), neighborsFromHorizontalEdge[k])];
      }
      if (update==Replace) {
        dst[indexFromHorizontalEdge<Level>(it.col(), it.row(), stencilDirection::EDGE_HO_C)] = tmp;
      } else if ( update==Add ) {
        dst[indexFromHorizontalEdge<Level>(it.col(), it.row(), stencilDirection::EDGE_HO_C)] += tmp;
      }
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k){
        tmp += opr_data[edgedof::stencilIndexFromDiagonalEdge(neighborsFromDiagonalEdge[k])] *
               src[indexFromDiagonalEdge< Level >(it.col(), it.row(), neighborsFromDiagonalEdge[k])];
      }
      if (update==Replace) {
        dst[indexFromDiagonalEdge< Level >(it.col(), it.row(), stencilDirection::EDGE_DI_C)] = tmp;
      } else if ( update==Add ) {
        dst[indexFromDiagonalEdge<Level>(it.col(), it.row(), stencilDirection::EDGE_DI_C)] += tmp;
      }
    }
    if( it.col() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k){
        tmp += opr_data[edgedof::stencilIndexFromVerticalEdge(neighborsFromVerticalEdge[k])] *
               src[indexFromVerticalEdge< Level >(it.col(), it.row(), neighborsFromVerticalEdge[k])];
      }

      if (update==Replace) {
        dst[indexFromVerticalEdge< Level >(it.col(), it.row(), stencilDirection::EDGE_VE_C)] = tmp;
      } else if ( update==Add ) {
        dst[indexFromVerticalEdge<Level>(it.col(), it.row(), stencilDirection::EDGE_VE_C)] += tmp;
      }
    }
  }
}

SPECIALIZE(void, applyTmpl, apply)

}
}
}


