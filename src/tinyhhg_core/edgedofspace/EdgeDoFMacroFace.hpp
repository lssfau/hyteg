
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
  const Point3D diagonalMicroEdgeOffset   = horizontalMicroEdgeOffset + verticalMicroEdgeOffset;

  // For each edge dof type, only one side of the face must not be updated.
  // Therefore we:
  //   1. iterate over the inner face
  //   2. iterate over the two missing borders for each edge dof type (skipping the outer most dofs)
  //   3. updating the missing dof on the opposite of the edge that must not be updated

  for ( const auto & it : indexing::edgedof::macroface::Iterator( Level, 1 ) )
  {
    const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( ( it.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( it.row() * 2     ) * verticalMicroEdgeOffset );
    const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( ( it.col() * 2     ) * horizontalMicroEdgeOffset + ( it.row() * 2 + 1 ) * verticalMicroEdgeOffset );
    const Point3D diagonalMicroEdgePosition   = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

    for (uint_t k = 0; k < srcPtr.size(); ++k) {
      srcVectorHorizontal[k] = srcPtr[k][indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() )];
      srcVectorVertical[k] = srcPtr[k][indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() )];
      srcVectorDiagonal[k] = srcPtr[k][indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() )];
    }

    faceData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = expr( horizontalMicroEdgePosition, srcVectorHorizontal );
    faceData[ indexing::edgedof::macroface::verticalIndex< Level >  ( it.col(), it.row() ) ] = expr( verticalMicroEdgePosition, srcVectorVertical );
    faceData[ indexing::edgedof::macroface::diagonalIndex< Level >  ( it.col(), it.row() ) ] = expr( diagonalMicroEdgePosition, srcVectorDiagonal );
  }

  // At the bottom face border, only the horizontal edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT, 0, 1 ) )
  {
    const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( ( it.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( it.row() * 2     ) * verticalMicroEdgeOffset );
    const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( ( it.col() * 2     ) * horizontalMicroEdgeOffset + ( it.row() * 2 + 1 ) * verticalMicroEdgeOffset );
    const Point3D diagonalMicroEdgePosition   = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

    for (uint_t k = 0; k < srcPtr.size(); ++k) {
      srcVectorHorizontal[k] = srcPtr[k][indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() )];
      srcVectorVertical[k] = srcPtr[k][indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() )];
      srcVectorDiagonal[k] = srcPtr[k][indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() )];
    }

    faceData[ indexing::edgedof::macroface::verticalIndex< Level >  ( it.col(), it.row() ) ] = expr( verticalMicroEdgePosition, srcVectorVertical );
    faceData[ indexing::edgedof::macroface::diagonalIndex< Level >  ( it.col(), it.row() ) ] = expr( diagonalMicroEdgePosition, srcVectorDiagonal );
  }

  // At the left face border, only the vertical edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::LEFT_BOTTOM_TO_TOP, 0, 1 ) )
  {
    const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( ( it.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( it.row() * 2     ) * verticalMicroEdgeOffset );
    const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( ( it.col() * 2     ) * horizontalMicroEdgeOffset + ( it.row() * 2 + 1 ) * verticalMicroEdgeOffset );
    const Point3D diagonalMicroEdgePosition   = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

    for (uint_t k = 0; k < srcPtr.size(); ++k) {
      srcVectorHorizontal[k] = srcPtr[k][indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() )];
      srcVectorVertical[k] = srcPtr[k][indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() )];
      srcVectorDiagonal[k] = srcPtr[k][indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() )];
    }

    faceData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = expr( horizontalMicroEdgePosition, srcVectorHorizontal );
    faceData[ indexing::edgedof::macroface::diagonalIndex< Level >  ( it.col(), it.row() ) ] = expr( diagonalMicroEdgePosition, srcVectorDiagonal );
  }

  // At the diagonal face border, only the diagonal edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 0, 1 ) )
  {
    const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( ( it.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( it.row() * 2     ) * verticalMicroEdgeOffset );
    const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( ( it.col() * 2     ) * horizontalMicroEdgeOffset + ( it.row() * 2 + 1 ) * verticalMicroEdgeOffset );
    const Point3D diagonalMicroEdgePosition   = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

    faceData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = expr( horizontalMicroEdgePosition, srcVectorHorizontal );
    faceData[ indexing::edgedof::macroface::verticalIndex< Level >  ( it.col(), it.row() ) ] = expr( verticalMicroEdgePosition, srcVectorVertical );
  }

  // Missing DoF at the opposite of the edge that must not be updated
  const auto topLeftCorner     = indexing::edgedof::macroface::getTopLeftCorner< Level >();
  const auto bottomLeftCorner  = indexing::edgedof::macroface::getBottomLeftCorner< Level >();
  const auto bottomRightCorner = indexing::edgedof::macroface::getBottomRightCorner< Level >();

  const uint_t horizontalDoFTopLeftCornerIdx   = indexing::edgedof::macroface::horizontalIndex< Level >( topLeftCorner.col(), topLeftCorner.row() );
  const uint_t diagonalDoFBottomLeftCornerIdx  = indexing::edgedof::macroface::diagonalIndex< Level >( bottomLeftCorner.col(), bottomLeftCorner.row() );
  const uint_t verticalDoFBottomRightCornerIdx = indexing::edgedof::macroface::verticalIndex< Level >( bottomRightCorner.col(), bottomRightCorner.row() );

  const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( ( topLeftCorner.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( topLeftCorner.row() * 2     ) * verticalMicroEdgeOffset );
  const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( ( bottomRightCorner.col() * 2     ) * horizontalMicroEdgeOffset + ( bottomRightCorner.row() * 2 + 1 ) * verticalMicroEdgeOffset );
  const Point3D diagonalMicroEdgePosition   =   ( faceBottomLeftCoords + ( ( bottomLeftCorner.col() * 2 + 1 ) * horizontalMicroEdgeOffset + ( bottomLeftCorner.row() * 2     ) * verticalMicroEdgeOffset ) )
                                              + verticalMicroEdgeOffset;

  for (uint_t k = 0; k < srcPtr.size(); ++k) {
    srcVectorHorizontal[k] = srcPtr[k][ horizontalDoFTopLeftCornerIdx ];
    srcVectorVertical[k] = srcPtr[k][ verticalDoFBottomRightCornerIdx ];
    srcVectorDiagonal[k] = srcPtr[k][ diagonalDoFBottomLeftCornerIdx ];
  }

  faceData[ horizontalDoFTopLeftCornerIdx ] = expr( horizontalMicroEdgePosition, srcVectorHorizontal );
  faceData[ verticalDoFBottomRightCornerIdx ] = expr( verticalMicroEdgePosition, srcVectorVertical );
  faceData[ diagonalDoFBottomLeftCornerIdx ] = expr( diagonalMicroEdgePosition, srcVectorDiagonal );
}

SPECIALIZE_WITH_VALUETYPE( void, interpolateTmpl, interpolate );


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

    // At the bottom face border, only the horizontal edge dofs must not be updated
    for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT, 0, 1 ) )
    {
      dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ];
    }

    // At the left face border, only the vertical edge dofs must not be updated
    for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::LEFT_BOTTOM_TO_TOP, 0, 1 ) )
    {
      dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] += scalar * srcData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ];
    }

    // At the diagonal face border, only the diagonal edge dofs must not be updated
    for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 0, 1 ) )
    {
      dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] += scalar * srcData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ];
      dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   += scalar * srcData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ];
    }

    // Missing DoF at the opposite of the edge that must not be updated
    const auto topLeftCorner     = indexing::edgedof::macroface::getTopLeftCorner< Level >();
    const auto bottomLeftCorner  = indexing::edgedof::macroface::getBottomLeftCorner< Level >();
    const auto bottomRightCorner = indexing::edgedof::macroface::getBottomRightCorner< Level >();

    const uint_t horizontalDoFTopLeftCornerIdx   = indexing::edgedof::macroface::horizontalIndex< Level >( topLeftCorner.col(), topLeftCorner.row() );
    const uint_t diagonalDoFBottomLeftCornerIdx  = indexing::edgedof::macroface::diagonalIndex< Level >( bottomLeftCorner.col(), bottomLeftCorner.row() );
    const uint_t verticalDoFBottomRightCornerIdx = indexing::edgedof::macroface::verticalIndex< Level >( bottomRightCorner.col(), bottomRightCorner.row() );

    dstData[ horizontalDoFTopLeftCornerIdx ]   += scalar * srcData[ horizontalDoFTopLeftCornerIdx ];
    dstData[ verticalDoFBottomRightCornerIdx ] += scalar * srcData[ verticalDoFBottomRightCornerIdx ];
    dstData[ diagonalDoFBottomLeftCornerIdx ]  += scalar * srcData[ diagonalDoFBottomLeftCornerIdx ];

  }
}

SPECIALIZE_WITH_VALUETYPE( void, addTmpl, add );

template< typename ValueType, uint_t Level >
inline void assignTmpl( Face & face, const std::vector< ValueType > & scalars,
                        const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcIds,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );

  auto dstData = face.getData( dstId )->getPointer( Level );

  for ( const auto & it : indexing::edgedof::macroface::Iterator( Level, 1 ) )
  {
    dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = static_cast< ValueType >( 0 );
    dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   = static_cast< ValueType >( 0 );
    dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   = static_cast< ValueType >( 0 );
  }

  // At the bottom face border, only the horizontal edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT, 0, 1 ) )
  {
    dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   = static_cast< ValueType >( 0 );
    dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   = static_cast< ValueType >( 0 );
  }

  // At the left face border, only the vertical edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::LEFT_BOTTOM_TO_TOP, 0, 1 ) )
  {
    dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = static_cast< ValueType >( 0 );
    dstData[ indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() ) ]   = static_cast< ValueType >( 0 );
  }

  // At the diagonal face border, only the diagonal edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 0, 1 ) )
  {
    dstData[ indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() ) ] = static_cast< ValueType >( 0 );
    dstData[ indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() ) ]   = static_cast< ValueType >( 0 );
  }

  // Missing DoF at the opposite of the edge that must not be updated
  const auto topLeftCorner     = indexing::edgedof::macroface::getTopLeftCorner< Level >();
  const auto bottomLeftCorner  = indexing::edgedof::macroface::getBottomLeftCorner< Level >();
  const auto bottomRightCorner = indexing::edgedof::macroface::getBottomRightCorner< Level >();

  const uint_t horizontalDoFTopLeftCornerIdx   = indexing::edgedof::macroface::horizontalIndex< Level >( topLeftCorner.col(), topLeftCorner.row() );
  const uint_t diagonalDoFBottomLeftCornerIdx  = indexing::edgedof::macroface::diagonalIndex< Level >( bottomLeftCorner.col(), bottomLeftCorner.row() );
  const uint_t verticalDoFBottomRightCornerIdx = indexing::edgedof::macroface::verticalIndex< Level >( bottomRightCorner.col(), bottomRightCorner.row() );

  dstData[ horizontalDoFTopLeftCornerIdx ]   = static_cast< ValueType >( 0 );
  dstData[ verticalDoFBottomRightCornerIdx ] = static_cast< ValueType >( 0 );
  dstData[ diagonalDoFBottomLeftCornerIdx ]  = static_cast< ValueType >( 0 );

  addTmpl< ValueType, Level >( face, scalars, srcIds, dstId );
}

SPECIALIZE_WITH_VALUETYPE( void, assignTmpl, assign );

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

  // At the bottom face border, only the horizontal edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::BOTTOM_LEFT_TO_RIGHT, 0, 1 ) )
  {
    const uint_t diagonalIdx   = indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() );
    const uint_t verticalIdx   = indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() );

    scalarProduct += lhsData[ diagonalIdx ]   * rhsData[ diagonalIdx ];
    scalarProduct += lhsData[ verticalIdx ]   * rhsData[ verticalIdx ];
  }

  // At the left face border, only the vertical edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::LEFT_BOTTOM_TO_TOP, 0, 1 ) )
  {
    const uint_t horizontalIdx = indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() );
    const uint_t diagonalIdx   = indexing::edgedof::macroface::diagonalIndex< Level >( it.col(), it.row() );

    scalarProduct += lhsData[ horizontalIdx ] * rhsData[ horizontalIdx ];
    scalarProduct += lhsData[ diagonalIdx ]   * rhsData[ diagonalIdx ];
  }

  // At the diagonal face border, only the diagonal edge dofs must not be updated
  for ( const auto & it : indexing::edgedof::macroface::BorderIterator( Level, indexing::FaceBorderDirection::DIAGONAL_BOTTOM_TO_TOP, 0, 1 ) )
  {
    const uint_t horizontalIdx = indexing::edgedof::macroface::horizontalIndex< Level >( it.col(), it.row() );
    const uint_t verticalIdx   = indexing::edgedof::macroface::verticalIndex< Level >( it.col(), it.row() );

    scalarProduct += lhsData[ horizontalIdx ] * rhsData[ horizontalIdx ];
    scalarProduct += lhsData[ verticalIdx ]   * rhsData[ verticalIdx ];
  }

  // Missing DoF at the opposite of the edge that must not be updated
  const auto topLeftCorner     = indexing::edgedof::macroface::getTopLeftCorner< Level >();
  const auto bottomLeftCorner  = indexing::edgedof::macroface::getBottomLeftCorner< Level >();
  const auto bottomRightCorner = indexing::edgedof::macroface::getBottomRightCorner< Level >();

  const uint_t horizontalDoFTopLeftCornerIdx   = indexing::edgedof::macroface::horizontalIndex< Level >( topLeftCorner.col(), topLeftCorner.row() );
  const uint_t diagonalDoFBottomLeftCornerIdx  = indexing::edgedof::macroface::diagonalIndex< Level >( bottomLeftCorner.col(), bottomLeftCorner.row() );
  const uint_t verticalDoFBottomRightCornerIdx = indexing::edgedof::macroface::verticalIndex< Level >( bottomRightCorner.col(), bottomRightCorner.row() );

  scalarProduct += lhsData[ horizontalDoFTopLeftCornerIdx ] * rhsData[ horizontalDoFTopLeftCornerIdx ];
  scalarProduct += lhsData[ diagonalDoFBottomLeftCornerIdx ] * rhsData[ diagonalDoFBottomLeftCornerIdx ];
  scalarProduct += lhsData[ verticalDoFBottomRightCornerIdx ] * rhsData[ verticalDoFBottomRightCornerIdx ];

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


