
#pragma once


#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>
#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;

template< typename ValueType >
inline void interpolate(const uint_t & Level, Face & face,
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

  Point3D xBlend;

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
        srcVectorHorizontal[k] = srcPtr[k][edgedof::macroface::horizontalIndex( Level, it.col(), it.row())];
      }

      face.getGeometryMap()->evalF(horizontalMicroEdgePosition, xBlend);
      faceData[edgedof::macroface::horizontalIndex( Level, it.col(), it.row())] = expr( xBlend, srcVectorHorizontal );
    }

    // Do not update vertical DoFs at left border
    if ( it.col() != 0 )
    {
      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
        srcVectorVertical[k] = srcPtr[k][edgedof::macroface::verticalIndex( Level, it.col(), it.row())];
      }

      face.getGeometryMap()->evalF(verticalMicroEdgePosition, xBlend);
      faceData[edgedof::macroface::verticalIndex( Level, it.col(), it.row())] = expr( xBlend, srcVectorVertical );
    }

    // Do not update diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      for ( uint_t k = 0; k < srcPtr.size(); ++k )
      {
        srcVectorDiagonal[k] = srcPtr[k][edgedof::macroface::diagonalIndex( Level, it.col(), it.row())];
      }

      face.getGeometryMap()->evalF(diagonalMicroEdgePosition, xBlend);
      faceData[edgedof::macroface::diagonalIndex( Level, it.col(), it.row())] = expr( xBlend, srcVectorDiagonal );
    }
  }
}



template< typename ValueType >
inline void add( const uint_t & Level, Face & face, const std::vector< ValueType > & scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );
  WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" );

  auto dstData = face.getData( dstId )->getPointer( Level );

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    ValueType tmpHorizontal = static_cast< ValueType >( 0.0 );
    ValueType tmpVertical   = static_cast< ValueType >( 0.0 );
    ValueType tmpDiagonal   = static_cast< ValueType >( 0.0 );

    const uint_t idxHorizontal = edgedof::macroface::horizontalIndex( Level, it.col(), it.row());
    const uint_t idxVertical   = edgedof::macroface::verticalIndex( Level, it.col(), it.row());
    const uint_t idxDiagonal   = edgedof::macroface::diagonalIndex( Level, it.col(), it.row());

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = face.getData( srcIds[i] )->getPointer( Level );

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
        tmpHorizontal += scalar * srcData[ idxHorizontal ];
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
        tmpVertical += scalar * srcData[ idxVertical ];
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
        tmpDiagonal += scalar * srcData[ idxDiagonal ];
      }
    }

    // Do not update horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      dstData[ idxHorizontal ] += tmpHorizontal;
    }

    // Do not update vertical DoFs at left border
    if ( it.col() != 0 )
    {
      dstData[ idxVertical ] += tmpVertical;
    }

    // Do not update diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      dstData[ idxDiagonal ] += tmpDiagonal;
    }
  }

}


template< typename ValueType >
inline void assign( const uint_t & Level, Face & face, const std::vector< ValueType > & scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > & srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstId )
{
  WALBERLA_ASSERT_EQUAL( scalars.size(), srcIds.size(), "Number of scalars must match number of src functions!" );
  WALBERLA_ASSERT_GREATER( scalars.size(), 0, "At least one src function and scalar must be given!" );

  auto dstData = face.getData( dstId )->getPointer( Level );

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    ValueType tmpHorizontal = static_cast< ValueType >( 0.0 );
    ValueType tmpVertical   = static_cast< ValueType >( 0.0 );
    ValueType tmpDiagonal   = static_cast< ValueType >( 0.0 );

    const uint_t idxHorizontal = edgedof::macroface::horizontalIndex( Level, it.col(), it.row());
    const uint_t idxVertical   = edgedof::macroface::verticalIndex( Level, it.col(), it.row());
    const uint_t idxDiagonal   = edgedof::macroface::diagonalIndex( Level, it.col(), it.row());

    for ( uint_t i = 0; i < scalars.size(); i++ )
    {
      const real_t scalar  = scalars[i];
      const auto   srcData = face.getData( srcIds[i] )->getPointer( Level );

      // Do not update horizontal DoFs at bottom
      if ( it.row() != 0 )
      {
        tmpHorizontal += scalar * srcData[ idxHorizontal ];
      }

      // Do not update vertical DoFs at left border
      if ( it.col() != 0 )
      {
        tmpVertical += scalar * srcData[ idxVertical ];
      }

      // Do not update diagonal DoFs at diagonal border
      if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
      {
        tmpDiagonal += scalar * srcData[ idxDiagonal ];
      }
    }

    // Do not update horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      dstData[ idxHorizontal ] = tmpHorizontal;
    }

    // Do not update vertical DoFs at left border
    if ( it.col() != 0 )
    {
      dstData[ idxVertical ] = tmpVertical;
    }

    // Do not update diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      dstData[ idxDiagonal ] = tmpDiagonal;
    }
  }
}


template< typename ValueType >
inline real_t dot( const uint_t & Level, Face & face,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& lhsId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId )
{
  auto lhsData = face.getData( lhsId )->getPointer( Level );
  auto rhsData = face.getData( rhsId )->getPointer( Level );

  walberla::math::KahanAccumulator< ValueType > scalarProduct;

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    // Do not read horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row());
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    // Do not read vertical DoFs at left border
    if ( it.col() != 0 )
    {
      const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row());
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }

    // Do not read diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row());
      scalarProduct += lhsData[ idx ] * rhsData[ idx ];
    }
  }

  return scalarProduct.get();
}


template< typename ValueType >
inline void enumerate(const uint_t & Level, Face &face,
                      const PrimitiveDataID < FunctionMemory< ValueType >, Face> &dstId,
                      uint_t& num)
{
  ValueType *dst = face.getData(dstId)->getPointer(Level);
  size_t horizontal_num = num;
  size_t diagonal_num = num +
                        hhg::edgedof::levelToFaceSizeAnyEdgeDoF( Level ) -
                        hhg::levelinfo::num_microedges_per_edge( Level ) ;
  size_t vertical_num = num +
                        (hhg::edgedof::levelToFaceSizeAnyEdgeDoF( Level ) -
                        hhg::levelinfo::num_microedges_per_edge( Level ))  *
                        2;
  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    /// the border edge DoFs belong to the corresponding edges
    if( it.row() != 0) {
      dst[hhg::edgedof::macroface::horizontalIndex( Level, it.col(), it.row())] = horizontal_num;
      ++horizontal_num;
      ++num;
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      dst[hhg::edgedof::macroface::diagonalIndex( Level, it.col(), it.row())] = diagonal_num;
      ++diagonal_num;
      ++num;
    }
    if( it.col() != 0) {
      dst[hhg::edgedof::macroface::verticalIndex( Level, it.col(), it.row())] = vertical_num;
      ++vertical_num;
      ++num;
    }

  }

}


inline void apply( const uint_t & Level, Face &face,
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
               src[indexFromHorizontalEdge( Level, it.col(), it.row(), neighborsFromHorizontalEdge[k] )];
      }
      if (update==Replace) {
        dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] = tmp;
      } else if ( update==Add ) {
        dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] += tmp;
      }
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k){
        tmp += opr_data[edgedof::stencilIndexFromDiagonalEdge(neighborsFromDiagonalEdge[k])] *
               src[indexFromDiagonalEdge( Level, it.col(), it.row(), neighborsFromDiagonalEdge[k] )];
      }
      if (update==Replace) {
        dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] = tmp;
      } else if ( update==Add ) {
        dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] += tmp;
      }
    }
    if( it.col() != 0) {
      tmp = 0.0;
      for(uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k){
        tmp += opr_data[edgedof::stencilIndexFromVerticalEdge(neighborsFromVerticalEdge[k])] *
               src[indexFromVerticalEdge( Level, it.col(), it.row(), neighborsFromVerticalEdge[k] )];
      }

      if (update==Replace) {
        dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] = tmp;
      } else if ( update==Add ) {
        dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] += tmp;
      }
    }
  }
}


template< typename ValueType >
inline void printFunctionMemory( const uint_t & Level, Face& face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId){
  ValueType* faceMemory = face.getData(dstId)->getPointer( Level );
  using namespace std;
  cout << setfill('=') << setw(100) << "" << endl;
  cout << face << std::left << setprecision(1) << fixed << setfill(' ') << endl;
  cout << "Horizontal Edge";
  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) ){
    if(it.col() == 0) std::cout << std::endl;
    cout << setw(5) << faceMemory[hhg::edgedof::macroface::indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )] << "|";
  }
  cout << endl << "Diagonal Edge";
  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) ){
    if(it.col() == 0) std::cout << std::endl;
    cout << setw(5) << faceMemory[hhg::edgedof::macroface::indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )] << "|";
  }
  cout << endl << "Vertical Edge";
  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) ){
    if(it.col() == 0) std::cout << std::endl;
    cout << setw(5) << faceMemory[hhg::edgedof::macroface::indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )] << "|";
  }
  cout << endl << setfill('=') << setw(100) << "" << endl << setfill(' ');

}


template< typename ValueType >
inline real_t getMaxMagnitude( const uint_t &level, Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId ) {

  auto src = face.getData( srcId )->getPointer( level );
  real_t localMax = real_t(0.0);

  for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
  {
    // Do not read horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      const uint_t idx = edgedof::macroface::horizontalIndex( level, it.col(), it.row() );
      localMax = std::max( localMax, std::abs( src[idx] ) );
    }

    // Do not read vertical DoFs at left border
    if ( it.col() != 0 )
    {
      const uint_t idx = edgedof::macroface::verticalIndex( level, it.col(), it.row() );
      localMax = std::max( localMax, std::abs( src[idx] ) );
    }

    // Do not read diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
    {
      const uint_t idx = edgedof::macroface::diagonalIndex( level, it.col(), it.row() );
      localMax = std::max( localMax, std::abs( src[idx] ) );
    }
  }

  return localMax;
}


#ifdef HHG_BUILD_WITH_PETSC

template< typename ValueType >
inline void createVectorFromFunction( const uint_t & Level, Face &face,
                                      const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                      const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &numeratorId,
                                      Vec& vec) {

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    // Do not read horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
      VecSetValues(vec,1,&numerator[idx],&src[idx],INSERT_VALUES);
    }

    // Do not read vertical DoFs at left border
    if ( it.col() != 0 )
    {
      const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
      VecSetValues(vec,1,&numerator[idx],&src[idx],INSERT_VALUES);
    }

    // Do not read diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );
      VecSetValues(vec,1,&numerator[idx],&src[idx],INSERT_VALUES);
    }
  }

}



template< typename ValueType >
inline void createFunctionFromVector( const uint_t & Level, Face &face,
                                      const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                      const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &numeratorId,
                                      Vec& vec) {

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );

  for ( const auto & it : edgedof::macroface::Iterator( Level, 0 ) )
  {
    // Do not read horizontal DoFs at bottom
    if ( it.row() != 0 )
    {
      const uint_t idx = edgedof::macroface::horizontalIndex( Level, it.col(), it.row() );
      VecGetValues(vec,1,&numerator[idx],&src[idx]);
    }

    // Do not read vertical DoFs at left border
    if ( it.col() != 0 )
    {
      const uint_t idx = edgedof::macroface::verticalIndex( Level, it.col(), it.row() );
      VecGetValues(vec,1,&numerator[idx],&src[idx]);
    }

    // Do not read diagonal DoFs at diagonal border
    if ( it.col() + it.row() != ( hhg::levelinfo::num_microedges_per_edge( Level ) - 1 ) )
    {
      const uint_t idx = edgedof::macroface::diagonalIndex( Level, it.col(), it.row() );
      VecGetValues(vec,1,&numerator[idx],&src[idx]);
    }
  }

}

#endif

}
}
}


