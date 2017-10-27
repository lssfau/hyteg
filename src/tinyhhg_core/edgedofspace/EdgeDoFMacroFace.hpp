
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

  for ( const auto & it : indexing::edgedof::macroface::Iterator< Level, 1 >() )
  {
    const Point3D horizontalMicroEdgePosition = ( it.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( it.row() * 2     ) * verticalMicroEdgeOffset;
    const Point3D verticalMicroEdgePosition   = ( it.col() * 2     ) * horizontalMicroEdgeOffset + ( it.row() * 2 + 1 )* verticalMicroEdgeOffset;
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

    for ( const auto & it : indexing::edgedof::macroface::Iterator< Level, 1 >() )
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

    for ( const auto & it : indexing::edgedof::macroface::Iterator< Level, 1 >() )
    {
      dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] += scalar * srcData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ];
    }
  }

}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add );

}
}
}


