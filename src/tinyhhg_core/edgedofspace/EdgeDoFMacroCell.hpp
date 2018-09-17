
#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/primitives/Cell.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"

namespace hhg {
namespace edgedof {
namespace macrocell {

using walberla::uint_t;
using walberla::real_c;

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


}
}
}