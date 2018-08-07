
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

}
}
}