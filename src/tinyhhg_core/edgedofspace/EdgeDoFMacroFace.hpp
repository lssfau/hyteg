
#pragma once


#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/macros.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMemory.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;

template< typename ValueType, uint_t Level >
inline void interpolateTmpl(Face & face,
                            const PrimitiveDataID< FunctionMemory< ValueType >, Face > & faceMemoryId,
                            std::function< ValueType( const hhg::Point3D & ) > & expr)
{
  auto faceData = face.getData( faceMemoryId )->getPointer( Level );

  const Point3D faceBottomLeftCoords  = face.coords[0];
  const Point3D faceBottomRightCoords = face.coords[1];
  const Point3D faceTopLeftCoords     = face.coords[2];

  const Point3D horizontalMicroEdgeOffset = ( ( faceBottomRightCoords - faceBottomLeftCoords ) / levelinfo::num_microedges_per_edge( Level ) ) * 0.5;
  const Point3D verticalMicroEdgeOffset   = ( ( faceTopLeftCoords     - faceBottomLeftCoords ) / levelinfo::num_microedges_per_edge( Level ) ) * 0.5;
  const Point3D diagonalMicroEdgeOffset   = horizontalMicroEdgeOffset + verticalMicroEdgeOffset;

  for ( const auto & it : indexing::edgedof::macroface::Iterator( Level, 1 ) )
  {
    const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( ( it.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( it.row() * 2     ) * verticalMicroEdgeOffset );
    const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( ( it.col() * 2     ) * horizontalMicroEdgeOffset + ( it.row() * 2 + 1 ) * verticalMicroEdgeOffset );
    const Point3D diagonalMicroEdgePosition   = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

    faceData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = expr( horizontalMicroEdgePosition );
    faceData[ indexing::edgedof::macroface::verticalIndex< Level >  ( it.col(), it.row() ) ] = expr( verticalMicroEdgePosition );
    faceData[ indexing::edgedof::macroface::diagonalIndex< Level >  ( it.col(), it.row() ) ] = expr( diagonalMicroEdgePosition );
  }
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate );

template< typename ValueType, uint_t Level >
inline void assignTmpl( Face & face, const std::vector< ValueType > & scalars,
                        const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcIds,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );

  auto dstData = face.getData( dstId )->getPointer( Level );

  for ( uint_t i = 0; i < scalars.size(); i++ )
  {
    const real_t scalar  = scalars[i];
    auto         srcData = face.getData( srcIds[i] )->getPointer( Level );

    for ( const auto & it : indexing::edgedof::macroface::Iterator( Level, 1 ) )
    {
      dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ]  = scalar * srcData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ];
    }
  }

}

SPECIALIZE_WITH_VALUETYPE( void, assignTmpl, assign );

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

    for ( const auto & it : indexing::edgedof::macroface::Iterator( Level, 1 ) )
    {
      dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] += scalar * srcData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ];
    }
  }

}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add );

template< typename ValueType, uint_t Level >
inline real_t dotTmpl( Face & face,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& lhsId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId )
{
  auto lhsData = face.getData( lhsId )->getPointer( Level );
  auto rhsData = face.getData( rhsId )->getPointer( Level );

  real_t scalarProduct = real_c( 0 );

  for ( const auto & it : indexing::edgedof::macroface::Iterator( Level, 1 ) )
  {
    const uint_t horizontalIdx = indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() );
    const uint_t diagonalIdx   = indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() );
    const uint_t verticalIdx   = indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() );

    scalarProduct += lhsData[ horizontalIdx ] * rhsData[ horizontalIdx ];
    scalarProduct += lhsData[ diagonalIdx ]   * rhsData[ diagonalIdx ];
    scalarProduct += lhsData[ verticalIdx ]   * rhsData[ verticalIdx ];
  }

  return scalarProduct;
}

SPECIALIZE_WITH_VALUETYPE( real_t, dotTmpl, dot );

template< typename ValueType, uint_t Level >
inline void enumerateTmpl(Face &face,
                          const PrimitiveDataID < FunctionMemory< ValueType >, Face> &dstId,
                          uint_t& num)
{
  ValueType *dst = face.getData(dstId)->getPointer(Level);
  size_t horizontal_num = num;
  size_t diagonal_num = num +
                        hhg::indexing::edgedof::levelToFaceSizeAnyEdgeDoF< Level > -
                        hhg::levelinfo::num_microedges_per_edge( Level ) ;
  size_t vertical_num = num +
                        (hhg::indexing::edgedof::levelToFaceSizeAnyEdgeDoF< Level > -
                        hhg::levelinfo::num_microedges_per_edge( Level ))  *
                        2;
  for ( const auto & it : hhg::indexing::edgedof::macroface::Iterator( Level, 0 ) )
  {
    /// the border edge DoFs belong to the corresponding edges
    if( it.row() != 0) {
      dst[hhg::indexing::edgedof::macroface::horizontalIndex< Level >(it.col(), it.row())] = horizontal_num;
      ++horizontal_num;
      ++num;
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      dst[hhg::indexing::edgedof::macroface::diagonalIndex< Level >(it.col(), it.row())] = diagonal_num;
      ++diagonal_num;
      ++num;
    }
    if( it.col() != 0) {
      dst[hhg::indexing::edgedof::macroface::verticalIndex< Level >(it.col(), it.row())] = vertical_num;
      ++vertical_num;
      ++num;
    }

  }

}

SPECIALIZE_WITH_VALUETYPE( void, enumerateTmpl, enumerate );

}
}
}


