
#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/primitives/Cell.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements3D.hpp"
#include "tinyhhg_core/Algorithms.hpp"

namespace hhg {
namespace edgedof {
namespace macrocell {

using walberla::uint_t;
using walberla::real_c;

inline indexing::Index getIndexInNeighboringMacroEdge( const indexing::Index  & edgeDoFIndexInMacroCell,
                                                       const Cell             & cell,
                                                       const uint_t           & neighborEdgeID,
                                                       const PrimitiveStorage & storage,
                                                       const uint_t           & level )
{
  const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 2, 4 >(
  { cell.getEdgeLocalVertexToCellLocalVertexMaps().at(neighborEdgeID).at(0),
    cell.getEdgeLocalVertexToCellLocalVertexMaps().at(neighborEdgeID).at(1) } );

  const auto indexInMacroEdge = indexing::basisConversion( edgeDoFIndexInMacroCell, {0, 1, 2, 3},
                                                           localVertexIDsAtCell, levelinfo::num_microedges_per_edge( level ) );
  return indexInMacroEdge;
}

inline indexing::Index getIndexInNeighboringMacroEdgeXYZ( const indexing::Index  & edgeDoFIndexInMacroCell,
                                                          const Cell             & cell,
                                                          const uint_t           & neighborEdgeID,
                                                          const PrimitiveStorage & storage,
                                                          const uint_t           & level )
{
  const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 2, 4 >(
  { cell.getEdgeLocalVertexToCellLocalVertexMaps().at(neighborEdgeID).at(0),
    cell.getEdgeLocalVertexToCellLocalVertexMaps().at(neighborEdgeID).at(1) } );

  const auto indexInMacroEdge = indexing::basisConversion( edgeDoFIndexInMacroCell, {0, 1, 2, 3},
                                                           localVertexIDsAtCell, levelinfo::num_microedges_per_edge( level ) - 1 );
  return indexInMacroEdge;
}

inline edgedof::EdgeDoFOrientation getOrientationInNeighboringMacroEdge( const EdgeDoFOrientation & orientationInMacroCell,
                                                                          const Cell               & cell,
                                                                          const uint_t             & neighborEdgeID,
                                                                          const PrimitiveStorage   & storage )
{
  const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 2, 4 >(
  { cell.getEdgeLocalVertexToCellLocalVertexMaps().at(neighborEdgeID).at(0),
    cell.getEdgeLocalVertexToCellLocalVertexMaps().at(neighborEdgeID).at(1) } );

  const auto orientationInMacroEdge = edgedof::convertEdgeDoFOrientationCellToFace( orientationInMacroCell,
                                                                                    localVertexIDsAtCell.at(0),
                                                                                    localVertexIDsAtCell.at(1),
                                                                                    localVertexIDsAtCell.at(2) );
  return orientationInMacroEdge;
}


inline indexing::Index getIndexInNeighboringMacroFace( const indexing::Index  & edgeDoFIndexInMacroCell,
                                                       const Cell             & cell,
                                                       const uint_t           & neighborFaceID,
                                                       const PrimitiveStorage & storage,
                                                       const uint_t           & level )
{
  const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 3, 4 >(
    { cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(0),
      cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(1),
      cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(2) } );

  const auto indexInMacroFace = indexing::basisConversion( edgeDoFIndexInMacroCell, {0, 1, 2, 3},
                                                           localVertexIDsAtCell, levelinfo::num_microedges_per_edge( level ) );
  return indexInMacroFace;
}

inline indexing::Index getIndexInNeighboringMacroFaceXYZ( const indexing::Index  & edgeDoFIndexInMacroCell,
                                                          const Cell             & cell,
                                                          const uint_t           & neighborFaceID,
                                                          const PrimitiveStorage & storage,
                                                          const uint_t           & level )
{
  const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 3, 4 >(
  { cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(0),
    cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(1),
    cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(2) } );

  const auto indexInMacroFace = indexing::basisConversion( edgeDoFIndexInMacroCell, {0, 1, 2, 3},
                                                           localVertexIDsAtCell, levelinfo::num_microedges_per_edge( level ) - 1 );
  return indexInMacroFace;
}

inline edgedof::EdgeDoFOrientation getOrientattionInNeighboringMacroFace( const EdgeDoFOrientation & orientationInMacroCell,
                                                                          const Cell               & cell,
                                                                          const uint_t             & neighborFaceID,
                                                                          const PrimitiveStorage   & storage )
{
  const auto orientationInMacroFace = edgedof::convertEdgeDoFOrientationCellToFace( orientationInMacroCell,
                                                                                    cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(0),
                                                                                    cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(1),
                                                                                    cell.getFaceLocalVertexToCellLocalVertexMaps().at(neighborFaceID).at(2) );
  return orientationInMacroFace;
}

inline Point3D xShiftFromVertex( const uint_t & level, const Cell & cell )
{
  const real_t stepFrequency = 1.0 / real_c( levelinfo::num_microedges_per_edge( level ) );
  return 0.5 * ( cell.getCoordinates()[1] - cell.getCoordinates()[0] ) * stepFrequency;
}

inline Point3D yShiftFromVertex( const uint_t & level, const Cell & cell )
{
  const real_t  stepFrequency = 1.0 / real_c( levelinfo::num_microedges_per_edge( level ) );
  return 0.5 * ( cell.getCoordinates()[2] - cell.getCoordinates()[0] ) * stepFrequency;
}

inline Point3D zShiftFromVertex( const uint_t & level, const Cell & cell )
{
  const real_t  stepFrequency = 1.0 / real_c( levelinfo::num_microedges_per_edge( level ) );
  return 0.5 * ( cell.getCoordinates()[3] - cell.getCoordinates()[0] ) * stepFrequency;
}


template< typename ValueType >
inline void interpolate(const uint_t & Level, Cell & cell,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & cellMemoryId,
                        const ValueType & constant )
{

  auto cellData = cell.getData( cellMemoryId )->getPointer( Level );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    cellData[edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z())] = constant;
    cellData[edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z())] = constant;
    cellData[edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z())] = constant;
    cellData[edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z())] = constant;
    cellData[edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z())] = constant;
    cellData[edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z())] = constant;
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    cellData[edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z())] = constant;
  }
}


template< typename ValueType >
inline void interpolate(const uint_t & Level, Cell & cell,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & cellMemoryId,
                        const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Cell>> &srcIds,
                        const std::function< ValueType( const hhg::Point3D &, const std::vector<ValueType>& ) > & expr)
{

  auto cellData = cell.getData( cellMemoryId )->getPointer( Level );

  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(cell.getData(src)->getPointer( Level ));
  }

  std::vector<ValueType> srcVectorX(srcIds.size());
  std::vector<ValueType> srcVectorY(srcIds.size());
  std::vector<ValueType> srcVectorZ(srcIds.size());
  std::vector<ValueType> srcVectorXY(srcIds.size());
  std::vector<ValueType> srcVectorXZ(srcIds.size());
  std::vector<ValueType> srcVectorYZ(srcIds.size());
  std::vector<ValueType> srcVectorXYZ(srcIds.size());

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    const Point3D microVertexPosition = vertexdof::macrocell::coordinateFromIndex( Level, cell, it );
    const Point3D xMicroEdgePosition  = microVertexPosition + xShiftFromVertex( Level, cell );
    const Point3D yMicroEdgePosition  = microVertexPosition + yShiftFromVertex( Level, cell );
    const Point3D zMicroEdgePosition  = microVertexPosition + zShiftFromVertex( Level, cell );
    const Point3D xyMicroEdgePosition = microVertexPosition + xShiftFromVertex( Level, cell ) + yShiftFromVertex( Level, cell );
    const Point3D xzMicroEdgePosition = microVertexPosition + xShiftFromVertex( Level, cell ) + zShiftFromVertex( Level, cell );
    const Point3D yzMicroEdgePosition = microVertexPosition + yShiftFromVertex( Level, cell ) + zShiftFromVertex( Level, cell );

    for ( uint_t k = 0; k < srcPtr.size(); ++k )
    {
      srcVectorX[k]  = srcPtr[k][edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z() )];
      srcVectorY[k]  = srcPtr[k][edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z() )];
      srcVectorZ[k]  = srcPtr[k][edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z() )];
      srcVectorXY[k] = srcPtr[k][edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z() )];
      srcVectorXZ[k] = srcPtr[k][edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z() )];
      srcVectorYZ[k] = srcPtr[k][edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z() )];
    }

    cellData[edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z())] = expr( xMicroEdgePosition, srcVectorX );
    cellData[edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z())] = expr( yMicroEdgePosition, srcVectorY );
    cellData[edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z())] = expr( zMicroEdgePosition, srcVectorZ );
    cellData[edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z())] = expr( xyMicroEdgePosition, srcVectorXY );
    cellData[edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z())] = expr( xzMicroEdgePosition, srcVectorXZ );
    cellData[edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z())] = expr( yzMicroEdgePosition, srcVectorYZ );
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    const Point3D microVertexPosition = vertexdof::macrocell::coordinateFromIndex( Level, cell, it );
    const Point3D xyzMicroEdgePosition  = microVertexPosition + xShiftFromVertex( Level, cell ) + yShiftFromVertex( Level, cell ) + zShiftFromVertex( Level, cell );

    for ( uint_t k = 0; k < srcPtr.size(); ++k )
    {
      srcVectorXYZ[k]  = srcPtr[k][edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z() )];
    }

    cellData[edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z())] = expr( xyzMicroEdgePosition, srcVectorXYZ );
  }
}


template< typename ValueType >
inline void assign(const uint_t & Level, Cell & cell, const std::vector< ValueType > & scalars,
                   const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );
  WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" );

  auto dstData = cell.getData( dstId )->getPointer( Level );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    ValueType tmpX  = static_cast< ValueType >( 0.0 );
    ValueType tmpY  = static_cast< ValueType >( 0.0 );
    ValueType tmpZ  = static_cast< ValueType >( 0.0 );
    ValueType tmpXY = static_cast< ValueType >( 0.0 );
    ValueType tmpXZ = static_cast< ValueType >( 0.0 );
    ValueType tmpYZ = static_cast< ValueType >( 0.0 );

    const uint_t idxX  = edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxY  = edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxZ  = edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxXY = edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxXZ = edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxYZ = edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z() );

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = cell.getData( srcIds[i] )->getPointer( Level );

      tmpX  += scalar * srcData[ idxX ];
      tmpY  += scalar * srcData[ idxY ];
      tmpZ  += scalar * srcData[ idxZ ];
      tmpXY += scalar * srcData[ idxXY ];
      tmpXZ += scalar * srcData[ idxXZ ];
      tmpYZ += scalar * srcData[ idxYZ ];
    }

    dstData[ idxX ] = tmpX;
    dstData[ idxY ] = tmpY;
    dstData[ idxZ ] = tmpZ;
    dstData[ idxXY ] = tmpXY;
    dstData[ idxXZ ] = tmpXZ;
    dstData[ idxYZ ] = tmpYZ;
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    ValueType tmpXYZ = static_cast< ValueType >( 0.0 );
    const uint_t idxXYZ = edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z() );

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = cell.getData( srcIds[i] )->getPointer( Level );

      tmpXYZ += scalar * srcData[ idxXYZ ];
    }

    dstData[ idxXYZ ] = tmpXYZ;
  }
}


template< typename ValueType >
inline void add(const uint_t & Level, Cell & cell, const std::vector< ValueType > & scalars,
                   const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > & srcIds,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );
  WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" );

  auto dstData = cell.getData( dstId )->getPointer( Level );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    ValueType tmpX  = static_cast< ValueType >( 0.0 );
    ValueType tmpY  = static_cast< ValueType >( 0.0 );
    ValueType tmpZ  = static_cast< ValueType >( 0.0 );
    ValueType tmpXY = static_cast< ValueType >( 0.0 );
    ValueType tmpXZ = static_cast< ValueType >( 0.0 );
    ValueType tmpYZ = static_cast< ValueType >( 0.0 );

    const uint_t idxX  = edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxY  = edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxZ  = edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxXY = edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxXZ = edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z() );
    const uint_t idxYZ = edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z() );

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = cell.getData( srcIds[i] )->getPointer( Level );

      tmpX  += scalar * srcData[ idxX ];
      tmpY  += scalar * srcData[ idxY ];
      tmpZ  += scalar * srcData[ idxZ ];
      tmpXY += scalar * srcData[ idxXY ];
      tmpXZ += scalar * srcData[ idxXZ ];
      tmpYZ += scalar * srcData[ idxYZ ];
    }

    dstData[ idxX ] += tmpX;
    dstData[ idxY ] += tmpY;
    dstData[ idxZ ] += tmpZ;
    dstData[ idxXY ] += tmpXY;
    dstData[ idxXZ ] += tmpXZ;
    dstData[ idxYZ ] += tmpYZ;
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    ValueType tmpXYZ = static_cast< ValueType >( 0.0 );
    const uint_t idxXYZ = edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z() );

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = cell.getData( srcIds[i] )->getPointer( Level );

      tmpXYZ += scalar * srcData[ idxXYZ ];
    }

    dstData[ idxXYZ ] += tmpXYZ;
  }
}


template< typename ValueType >
inline real_t dot( const uint_t & Level, Cell & cell,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& lhsId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& rhsId )
{
  auto lhsData = cell.getData( lhsId )->getPointer( Level );
  auto rhsData = cell.getData( rhsId )->getPointer( Level );

  walberla::math::KahanAccumulator< ValueType > scalarProduct;

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    if ( isInnerXEdgeDoF( Level, it ) )
    {
      const uint_t idx  = edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    if ( isInnerYEdgeDoF( Level, it ) )
    {
      const uint_t idx  = edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    if ( isInnerZEdgeDoF( Level, it ) )
    {
      const uint_t idx  = edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    if ( isInnerXYEdgeDoF( Level, it ) )
    {
      const uint_t idx  = edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    if ( isInnerXZEdgeDoF( Level, it ) )
    {
      const uint_t idx  = edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    if ( isInnerYZEdgeDoF( Level, it ) )
    {
      const uint_t idx  = edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z() );
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    const uint_t idx  = edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z() );
    scalarProduct += lhsData[ idx ] * rhsData[ idx ];
  }

  return scalarProduct.get();
}

inline void apply(const uint_t & Level, Cell &cell,
                  const PrimitiveDataID< LevelWiseMemory< StencilMap_T >, Cell> &operatorId,
                  const PrimitiveDataID< FunctionMemory< real_t >, Cell> &srcId,
                  const PrimitiveDataID< FunctionMemory< real_t >, Cell> &dstId,
                  UpdateType update)
{
  auto srcData = cell.getData( srcId )->getPointer( Level );
  auto dstData = cell.getData( dstId )->getPointer( Level );
  auto opr_data = cell.getData( operatorId )->getData( Level );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

    if ( isInnerXEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
    if ( isInnerYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
    if ( isInnerZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
    if ( isInnerXYEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
    if ( isInnerXZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
    if ( isInnerYZEdgeDoF( Level, it ) )
      innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

    for ( const auto & centerOrientation : innerOrientations )
    {
      real_t tmp = 0.0;
      for ( const auto & leafOrientation : edgedof::allEdgeDoFOrientations )
      {
        const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
        for ( const auto & neighbor : edgeDoFNeighbors )
        {
          const auto   srcIdx      = it + neighbor;
          const auto   srcArrayIdx = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
          const real_t stencilWeight = opr_data[centerOrientation][leafOrientation][neighbor];
          tmp += stencilWeight * srcData[srcArrayIdx];
        }
      }

      const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );

      if ( update == Replace )
      {
        dstData[dstArrayIdx] = tmp;
      }
      else if ( update == Add )
      {
        dstData[dstArrayIdx] += tmp;
      }
    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    real_t tmp = 0.0;
    const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;
    for ( const auto & leafOrientation : edgedof::allEdgeDoFOrientations )
    {
      const auto edgeDoFNeighbors = P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
      for ( const auto & neighbor : edgeDoFNeighbors )
      {
        const auto   srcIdx      = it + neighbor;
        const auto   srcArrayIdx = edgedof::macrocell::index( Level, srcIdx.x(), srcIdx.y(), srcIdx.z(), leafOrientation );
        const real_t stencilWeight = opr_data[centerOrientation][leafOrientation][neighbor];
        tmp += stencilWeight * srcData[srcArrayIdx];
      }
    }

    const auto dstArrayIdx = edgedof::macrocell::index( Level, it.x(), it.y(), it.z(), centerOrientation );

    if ( update == Replace )
    {
      dstData[dstArrayIdx] = tmp;
    }
    else if ( update == Add )
    {
      dstData[dstArrayIdx] += tmp;
    }
  }
}

template< typename ValueType >
inline ValueType getMaxValue( const uint_t &level, Cell &cell, const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId ) {

  auto src = cell.getData( srcId )->getPointer( level );
  auto localMax = -std::numeric_limits< ValueType >::max();

  for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) ) {
    localMax = std::max( localMax, src[ edgedof::macrocell::xIndex( level, it.x(), it.y(), it.z() ) ] );
    localMax = std::max( localMax, src[ edgedof::macrocell::yIndex( level, it.x(), it.y(), it.z() ) ] );
    localMax = std::max( localMax, src[ edgedof::macrocell::zIndex( level, it.x(), it.y(), it.z() ) ] );
    localMax = std::max( localMax, src[ edgedof::macrocell::xyIndex( level, it.x(), it.y(), it.z() ) ] );
    localMax = std::max( localMax, src[ edgedof::macrocell::xzIndex( level, it.x(), it.y(), it.z() ) ] );
    localMax = std::max( localMax, src[ edgedof::macrocell::yzIndex( level, it.x(), it.y(), it.z() ) ] );
  }

  for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) ) {
    localMax = std::max( localMax, src[ edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() ) ] );
  }

  return localMax;
}

template< typename ValueType >
inline ValueType getMinValue( const uint_t &level, Cell &cell, const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId ) {

  auto src = cell.getData( srcId )->getPointer( level );
  auto localMin = std::numeric_limits< ValueType >::max();

  for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) ) {
    localMin = std::min( localMin, src[ edgedof::macrocell::xIndex( level, it.x(), it.y(), it.z() ) ] );
    localMin = std::min( localMin, src[ edgedof::macrocell::yIndex( level, it.x(), it.y(), it.z() ) ] );
    localMin = std::min( localMin, src[ edgedof::macrocell::zIndex( level, it.x(), it.y(), it.z() ) ] );
    localMin = std::min( localMin, src[ edgedof::macrocell::xyIndex( level, it.x(), it.y(), it.z() ) ] );
    localMin = std::min( localMin, src[ edgedof::macrocell::xzIndex( level, it.x(), it.y(), it.z() ) ] );
    localMin = std::min( localMin, src[ edgedof::macrocell::yzIndex( level, it.x(), it.y(), it.z() ) ] );
  }

  for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) ) {
    localMin = std::min( localMin, src[ edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() ) ] );
  }

  return localMin;
}

template< typename ValueType >
inline ValueType getMaxMagnitude( const uint_t &level, Cell &cell, const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId ) {

  auto src = cell.getData( srcId )->getPointer( level );
  auto localMax = ValueType(0.0);

  for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) ) {
    localMax = std::max( localMax, std::abs( src[ edgedof::macrocell::xIndex( level, it.x(), it.y(), it.z() ) ] ) );
    localMax = std::max( localMax, std::abs( src[ edgedof::macrocell::yIndex( level, it.x(), it.y(), it.z() ) ] ) );
    localMax = std::max( localMax, std::abs( src[ edgedof::macrocell::zIndex( level, it.x(), it.y(), it.z() ) ] ) );
    localMax = std::max( localMax, std::abs( src[ edgedof::macrocell::xyIndex( level, it.x(), it.y(), it.z() ) ] ) );
    localMax = std::max( localMax, std::abs( src[ edgedof::macrocell::xzIndex( level, it.x(), it.y(), it.z() ) ] ) );
    localMax = std::max( localMax, std::abs( src[ edgedof::macrocell::yzIndex( level, it.x(), it.y(), it.z() ) ] ) );
  }

  for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) ) {
    localMax = std::max( localMax, std::abs( src[ edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() ) ] ) );
  }

  return localMax;
}


template< typename ValueType >
inline void enumerate(const uint_t & Level, Cell &cell,
                      const PrimitiveDataID < FunctionMemory< ValueType >, Cell > & dstId,
                      ValueType& num)
{
  ValueType *dst = cell.getData(dstId)->getPointer(Level);
  uint_t xNum   = uint_c( num );
  uint_t yNum   = xNum  + levelinfo::num_microvertices_per_cell_from_width( ( uint_c(1) << Level ) - 2 );
  uint_t zNum   = yNum  + levelinfo::num_microvertices_per_cell_from_width( ( uint_c(1) << Level ) - 2 );
  uint_t xyNum  = zNum  + levelinfo::num_microvertices_per_cell_from_width( ( uint_c(1) << Level ) - 2 );
  uint_t xzNum  = xyNum + levelinfo::num_microvertices_per_cell_from_width( ( uint_c(1) << Level ) - 2 );
  uint_t yzNum  = xzNum + levelinfo::num_microvertices_per_cell_from_width( ( uint_c(1) << Level ) - 2 );
  uint_t xyzNum = yzNum + levelinfo::num_microvertices_per_cell_from_width( ( uint_c(1) << Level ) - 2 );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    if ( isInnerXEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z() );
      dst[ idx ] = xNum;
      xNum++;
      num++;
    }

    if ( isInnerYEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z() );
      dst[ idx ] = yNum;
      yNum++;
      num++;    }

    if ( isInnerZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z() );
      dst[ idx ] = zNum;
      zNum++;
      num++;    }

    if ( isInnerXYEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z() );
      dst[ idx ] = xyNum;
      xyNum++;
      num++;    }

    if ( isInnerXZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z() );
      dst[ idx ] = xzNum;
      xzNum++;
      num++;    }

    if ( isInnerYZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z() );
      dst[ idx ] = yzNum;
      yzNum++;
      num++;    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    const uint_t idx = edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z() );
    dst[ idx ] = xyzNum;
    xyzNum++;
    num++;
  }

}


#ifdef HHG_BUILD_WITH_PETSC

template< typename ValueType >
inline void createVectorFromFunction( const uint_t & Level, Cell & cell,
                                      const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &srcId,
                                      const PrimitiveDataID<FunctionMemory< PetscInt >, Cell> &numeratorId,
                                      Vec& vec) {

  auto src = cell.getData(srcId)->getPointer( Level );
  auto numerator = cell.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    if ( isInnerXEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z() );
      VecSetValues( vec, 1, &numerator[idx], &src[idx], INSERT_VALUES );
    }

    if ( isInnerYEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z() );
      VecSetValues( vec, 1, &numerator[idx], &src[idx], INSERT_VALUES );
    }

    if ( isInnerZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z() );
      VecSetValues( vec, 1, &numerator[idx], &src[idx], INSERT_VALUES );
    }

    if ( isInnerXYEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z() );
      VecSetValues( vec, 1, &numerator[idx], &src[idx], INSERT_VALUES );
    }

    if ( isInnerXZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z() );
      VecSetValues( vec, 1, &numerator[idx], &src[idx], INSERT_VALUES );
    }

    if ( isInnerYZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z() );
      VecSetValues( vec, 1, &numerator[idx], &src[idx], INSERT_VALUES );
    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    const uint_t idx = edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z() );
    VecSetValues( vec, 1, &numerator[idx], &src[idx], INSERT_VALUES );
  }
}


template< typename ValueType >
inline void createFunctionFromVector( const uint_t & Level, Cell & cell,
                                      const PrimitiveDataID<FunctionMemory< ValueType >, Cell> &dstId,
                                      const PrimitiveDataID<FunctionMemory< PetscInt >, Cell> &numeratorId,
                                      Vec& vec) {
  
  auto dst = cell.getData(dstId)->getPointer( Level );
  auto numerator = cell.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : edgedof::macrocell::Iterator( Level, 0 ) )
  {
    if ( isInnerXEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xIndex( Level, it.x(), it.y(), it.z() );
      VecGetValues( vec, 1, &numerator[idx], &dst[idx] );
    }

    if ( isInnerYEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::yIndex( Level, it.x(), it.y(), it.z() );
      VecGetValues( vec, 1, &numerator[idx], &dst[idx] );
    }

    if ( isInnerZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::zIndex( Level, it.x(), it.y(), it.z() );
      VecGetValues( vec, 1, &numerator[idx], &dst[idx] );
    }

    if ( isInnerXYEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xyIndex( Level, it.x(), it.y(), it.z() );
      VecGetValues( vec, 1, &numerator[idx], &dst[idx] );
    }

    if ( isInnerXZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::xzIndex( Level, it.x(), it.y(), it.z() );
      VecGetValues( vec, 1, &numerator[idx], &dst[idx] );
    }

    if ( isInnerYZEdgeDoF( Level, it ) )
    {
      const uint_t idx = edgedof::macrocell::yzIndex( Level, it.x(), it.y(), it.z() );
      VecGetValues( vec, 1, &numerator[idx], &dst[idx] );
    }
  }

  for ( const auto & it : edgedof::macrocell::IteratorXYZ( Level, 0 ) )
  {
    const uint_t idx = edgedof::macrocell::xyzIndex( Level, it.x(), it.y(), it.z() );
    VecGetValues( vec, 1, &numerator[idx], &dst[idx] );
  }
}

#endif

}
}
}
