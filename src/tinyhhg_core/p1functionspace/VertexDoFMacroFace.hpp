#pragma once

#include <algorithm>
#include <cstdlib>

#include "core/debug/all.h"
#include "core/math/KahanSummation.h"

#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/facedofspace/FaceDoFIndexing.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include "tinyhhg_core/indexing/DistanceCoordinateSystem.hpp"

#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/primitives/Cell.hpp"
#include "tinyhhg_core/Algorithms.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {

using indexing::Index;
using walberla::real_c;
using walberla::uint_t;

inline indexing::Index getIndexInNeighboringMacroCell( const indexing::Index  & vertexDoFIndexInMacroFace,
                                                       const Face             & face,
                                                       const uint_t           & neighborCellID,
                                                       const PrimitiveStorage & storage,
                                                       const uint_t           & level )
{
   const Cell & neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
   const uint_t localFaceID = neighborCell.getLocalFaceID( face.getID() );

   const std::array< uint_t, 4 > localVertexIDsAtCell = algorithms::getMissingIntegersAscending< 3, 4 >(
   { neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(0),
     neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(1),
     neighborCell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(2) } );

   const auto indexInMacroCell = indexing::basisConversion( vertexDoFIndexInMacroFace, localVertexIDsAtCell,
                                                            {0, 1, 2, 3}, levelinfo::num_microvertices_per_edge( level ) );
   return indexInMacroCell;
}

inline Point3D coordinateFromIndex( const uint_t& Level, const Face& face, const Index& index )
{
   const real_t  stepFrequency = 1.0 / levelinfo::num_microedges_per_edge( Level );
   const Point3D xStep         = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) * stepFrequency;
   const Point3D yStep         = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) * stepFrequency;
   return face.getCoordinates()[0] + xStep * real_c( index.x() ) + yStep * real_c( index.y() );
}

template < typename ValueType >
inline ValueType assembleLocal( const uint_t&                            Level,
                                uint_t                                   i,
                                uint_t                                   j,
                                const Matrix3r&                          localMatrix,
                                double*                                  src,
                                double*                                  coeff,
                                const std::array< stencilDirection, 3 >& vertices,
                                const std::array< uint_t, 3 >&           idx )
{
   ValueType meanCoeff = 1.0 / 3.0 *
                         ( coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )] );

   ValueType tmp;
   tmp = localMatrix( idx[0], idx[0] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
         localMatrix( idx[0], idx[1] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
         localMatrix( idx[0], idx[2] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )];
   return meanCoeff * tmp;
}

template < typename ValueType >
inline void assembleLocalStencil( uint_t                                   Level,
                                  uint_t                                   i,
                                  uint_t                                   j,
                                  const Matrix3r&                          localMatrix,
                                  double*                                  opr_data,
                                  double*                                  coeff,
                                  const std::array< stencilDirection, 3 >& vertices,
                                  const std::array< uint_t, 3 >&           idx )
{
   ValueType meanCoeff = 1.0 / 3.0 *
                         ( coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
                           coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )] );

   opr_data[vertexdof::stencilIndexFromVertex( vertices[0] )] += meanCoeff * localMatrix( idx[0], idx[0] );
   opr_data[vertexdof::stencilIndexFromVertex( vertices[1] )] += meanCoeff * localMatrix( idx[0], idx[1] );
   opr_data[vertexdof::stencilIndexFromVertex( vertices[2] )] += meanCoeff * localMatrix( idx[0], idx[2] );
}

template < typename ValueType >
inline ValueType assembleLocalDG( const uint_t&                            Level,
                                  uint_t                                   i,
                                  uint_t                                   j,
                                  const Matrix3r&                          localMatrix,
                                  double*                                  src,
                                  const std::array< stencilDirection, 3 >& vertices,
                                  const std::array< uint_t, 3 >&           idx )
{
   ValueType tmp;
   tmp = localMatrix( idx[0], idx[0] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )] +
         localMatrix( idx[0], idx[1] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )] +
         localMatrix( idx[0], idx[2] ) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )];
   return tmp;
}

template < typename ValueType >
inline void getLocalElementDoFIndicesFromCoordinates( const uint_t & level, Face & face, const Point3D coordinates,
                                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face > & srcID,
                                                      Point2D& localCoordinates, Matrix2r& transform, Point3D& dofs)
{
  // Get the element local coordinates and DoFs from physical coordinates
  // The local DoFs are sorted in following order
  // 2
  // | \
  // |  \
  // 0 - 1

  // Transform absolute coordinates to macro element relative coordinates
  Matrix2r A;
  A(0,0) = (face.getCoordinates()[1] - face.getCoordinates()[0])[0];
  A(0,1) = (face.getCoordinates()[2] - face.getCoordinates()[0])[0];
  A(1,0) = (face.getCoordinates()[1] - face.getCoordinates()[0])[1];
  A(1,1) = (face.getCoordinates()[2] - face.getCoordinates()[0])[1];
  transform = A.adj();
  transform *= 1.0 / A.det();

  Point2D x({coordinates[0] - face.getCoordinates()[0][0], coordinates[1] - face.getCoordinates()[0][1]});

  Point2D xRelMacro = transform.mul(x);

  // Determine lower-left corner index of the quad where the evaluation point lies in
  uint_t rowsize = levelinfo::num_microvertices_per_edge( level );
  real_t hInv = walberla::real_c(rowsize - 1);
  real_t h = walberla::real_c(1.0 / hInv);
  Index index;
  index.x() = walberla::uint_c(std::floor(xRelMacro[0] * (rowsize-1)));
  index.y() = walberla::uint_c(std::floor(xRelMacro[1] * (rowsize-1)));

  if (index.x() == rowsize-1) {
    --index.x();
  }

  if (index.y() == rowsize-1) {
    --index.y();
  }

  localCoordinates[0] = xRelMacro[0] - index.x() * h;
  localCoordinates[1] = xRelMacro[1] - index.y() * h;
  localCoordinates *= hInv;

  auto srcData = face.getData( srcID )->getPointer( level );

  transform *= hInv;
  transform = transform.transpose();

  // Combine solution from linear basis functions depending on triangle orientation
  // Up triangle
  if (localCoordinates[0] + localCoordinates[1] <= 1.0) {
    dofs[0] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_C)];
    dofs[1] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_E)];
    dofs[2] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_N)];
  } else { // Down triangle
    index.y() += 1;
    localCoordinates[0] = 1.0 - localCoordinates[0];
    localCoordinates[1] = 1.0 - localCoordinates[1];
    transform *= -1.0;
    dofs[0] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_E)];
    dofs[1] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_C)];
    dofs[2] = srcData[vertexdof::macroface::indexFromVertex( level, index.x(), index.y(), stencilDirection::VERTEX_SE)];
  }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                               Level,
                         Face&                                                       face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >& faceMemoryId,
                         const ValueType&                                            scalar )
{
   ValueType* faceData = face.getData( faceMemoryId )->getPointer( Level );

   for( const auto& it : vertexdof::macroface::Iterator( Level, 1 ) )
   {
      const uint_t idx = vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(), stencilDirection::VERTEX_C );
      faceData[idx]    = scalar;
   }
}

template < typename ValueType >
inline void interpolate( const uint_t&                                                                             Level,
                         Face&                                                                                     face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                               faceMemoryId,
                         const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >&                srcIds,
                         const std::function< ValueType( const hhg::Point3D&, const std::vector< ValueType >& ) >& expr )
{
   ValueType* faceData = face.getData( faceMemoryId )->getPointer( Level );

   std::vector< ValueType* > srcPtr;
   for( const auto& src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   std::vector< ValueType > srcVector( srcIds.size() );

   Point3D xBlend;

   for( const auto& it : vertexdof::macroface::Iterator( Level, 1 ) )
   {
      const Point3D coordinate = coordinateFromIndex( Level, face, it );
      const uint_t  idx        = vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(), stencilDirection::VERTEX_C );

      for( uint_t k = 0; k < srcPtr.size(); ++k )
      {
         srcVector[k] = srcPtr[k][idx];
      }
      face.getGeometryMap()->evalF( coordinate, xBlend );
      faceData[idx] = expr( xBlend, srcVector );
   }
}

template< typename ValueType >
inline real_t evaluate( const uint_t & level, Face & face, const Point3D coordinates,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face > & srcID)
{
   Point2D localCoordinates;
   Matrix2r transform;
   Point3D localDoFs;
   getLocalElementDoFIndicesFromCoordinates( level, face, coordinates, srcID, localCoordinates, transform, localDoFs );
   real_t value = (1.0 - localCoordinates[0] - localCoordinates[1]) * localDoFs[0];
   value += localCoordinates[0] * localDoFs[1];
   value += localCoordinates[1] * localDoFs[2];

   return value;
}

template< typename ValueType >
inline void evaluateGradient( const uint_t & level, Face & face, const Point3D coordinates,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face > & srcID, Point3D& gradient)
{
  Point2D localCoordinates;
  Matrix2r transform;
  Point3D localDoFs;
  getLocalElementDoFIndicesFromCoordinates( level, face, coordinates, srcID, localCoordinates, transform, localDoFs );
  Point2D gradient_;
  gradient_[0] = localDoFs[1] - localDoFs[0];
  gradient_[1] = localDoFs[2] - localDoFs[0];
  gradient_ = transform.mul(gradient_);
  gradient[0] = gradient_[0];
  gradient[1] = gradient_[1];
}

template< typename ValueType >
inline void swap( const uint_t & level, Face & face,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & srcID,
                  const PrimitiveDataID< FunctionMemory< ValueType >, Face > & dstID )
{
  auto srcData = face.getData( srcID );
  auto dstData = face.getData( dstID );
  srcData->swap( *dstData, level );
}


template < typename ValueType >
inline void assign( const uint_t&                                                              Level,
                    Face&                                                                      face,
                    const std::vector< ValueType >&                                            scalars,
                    const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType*                dst = face.getData( dstId )->getPointer( Level );
   std::vector< ValueType* > srcPtr;
   for( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }
   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         ValueType tmp = scalars[0] * srcPtr[0][vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         for( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * srcPtr[k][vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         }
         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                               level,
                 const Face&                                                 face,
                 const ValueType&                                            scalar,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   ValueType* dstPtr = face.getData( dstId )->getPointer( level );

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         dstPtr[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] += scalar;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void add( const uint_t&                                                              Level,
                 Face&                                                                      face,
                 const std::vector< ValueType >&                                            scalars,
                 const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                 const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType*                dstPtr = face.getData( dstId )->getPointer( Level );
   std::vector< ValueType* > srcPtr;
   for( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         ValueType tmp = 0.0;

         for( uint_t k = 0; k < srcIds.size(); ++k )
         {
            tmp += scalars[k] * srcPtr[k][vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         }
         dstPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
      }
      --inner_rowsize;
   }
}


template < typename ValueType >
inline void multElementwise( const uint_t&                                                              Level,
                             Face&                                                                      face,
                             const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcIds,
                             const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType*                dstPtr = face.getData( dstId )->getPointer( Level );
   std::vector< ValueType* > srcPtr;
   for( auto src : srcIds )
   {
      srcPtr.push_back( face.getData( src )->getPointer( Level ) );
   }

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         const uint_t idx = vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C );
         ValueType tmp = srcPtr[0][idx];

         for( uint_t k = 1; k < srcIds.size(); ++k )
         {
            tmp *= srcPtr[k][idx];
         }
         dstPtr[idx] = tmp;
      }
      --inner_rowsize;
   }
}


template < typename ValueType >
inline real_t dot( const uint_t&                                               Level,
                   Face&                                                       face,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& lhsId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId )
{
   walberla::math::KahanAccumulator< ValueType > scalarProduct;
   uint_t                                        rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t                                        inner_rowsize = rowsize;

   ValueType* lhsPtr = face.getData( lhsId )->getPointer( Level );
   ValueType* rhsPtr = face.getData( rhsId )->getPointer( Level );

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         scalarProduct += lhsPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] *
                          rhsPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
      }
      --inner_rowsize;
   }

   return scalarProduct.get();
}

template < typename ValueType >
inline real_t sum( const uint_t& level, const Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dataID )
{
   ValueType* faceData = face.getData( dataID )->getPointer( level );
   real_t     sum      = real_c( 0 );
   for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
   {
      const uint_t idx = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C );
      sum += faceData[idx];
   }
   return sum;
}

template < typename ValueType >
inline void apply( const uint_t&                                               Level,
                   Face&                                                       face,
                   const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                   UpdateType                                                  update )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType* opr_data = face.getData( operatorId )->getPointer( Level );
   ValueType* src      = face.getData( srcId )->getPointer( Level );
   ValueType* dst      = face.getData( dstId )->getPointer( Level );

   ValueType tmp = real_c( 0 );

  for( uint_t j = 1; j < rowsize - 2; ++j )
  {
     for( uint_t i = 1; i < inner_rowsize - 2; ++i )
     {
        if( face.getNumNeighborCells() == 0 )
        {
           tmp = opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] *
                 src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

           //strangely the intel compiler cant handle this if it is a loop
           static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
           tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[0] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0] )];
           tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[1] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1] )];
           tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[2] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2] )];
           tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[3] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3] )];
           tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[4] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4] )];
           tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[5] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5] )];
        } else if( face.getNumNeighborCells() == 1 )
        {
           tmp = real_c( 0 );
           for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter )
           {
              tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                     src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
           }
        } else if( face.getNumNeighborCells() == 2 )
        {
           tmp = real_c( 0 );
           for( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter )
           {
              tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                     src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
           }
        }

        WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

        if ( update == Replace )
        {
          dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
        }
        else
        {
          dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
        }
     }
     --inner_rowsize;
  }
}


template < typename ValueType >
inline void applyPointwise( const uint_t&                                               Level,
                   Face&                                                       face,
                   const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                   UpdateType                                                  update )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   ValueType* opr_data = face.getData( operatorId )->getPointer( Level );
   ValueType* src      = face.getData( srcId )->getPointer( Level );
   ValueType* dst      = face.getData( dstId )->getPointer( Level );

   ValueType tmp = real_c( 0 );

  const uint_t stencilSize = 7;

  for( uint_t j = 1; j < rowsize - 2; ++j )
  {
     for( uint_t i = 1; i < inner_rowsize - 2; ++i )
     {
        if( face.getNumNeighborCells() == 0 )
        {
          const auto offset = stencilSize * vertexdof::macroface::innerIndex( Level, i, j );

           tmp = opr_data[offset + vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] *
                 src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

           //strangely the intel compiler cant handle this if it is a loop
           static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
           tmp += opr_data[offset + vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[0] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0] )];
           tmp += opr_data[offset + vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[1] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1] )];
           tmp += opr_data[offset + vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[2] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2] )];
           tmp += opr_data[offset + vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[3] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3] )];
           tmp += opr_data[offset + vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[4] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4] )];
           tmp += opr_data[offset + vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[5] )] *
                  src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5] )];
        } else if( face.getNumNeighborCells() == 1 )
        {
           tmp = real_c( 0 );
           for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter )
           {
              tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                     src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
           }
        } else if( face.getNumNeighborCells() == 2 )
        {
           tmp = real_c( 0 );
           for( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter )
           {
              tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                     src[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
           }
        }

        WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

        if ( update == Replace )
        {
          dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
        }
        else
        {
          dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
        }
     }
     --inner_rowsize;
  }
}

template < typename ValueType >
inline void applyCoefficient( const uint_t&                                                              Level,
                              Face&                                                                      face,
                              const std::vector< PrimitiveDataID< FaceP1LocalMatrixMemory, Face > >&     operatorIds,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                srcId,
                              const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId,
                              const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& coeffIds,
                              UpdateType                                                                 update )
{
   typedef stencilDirection SD;

   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto src = face.getData( srcId )->getPointer( Level );
   auto dst = face.getData( dstId )->getPointer( Level );

   std::vector< FaceP1LocalMatrixMemory* > localMatricesVector;
   for( auto operatorId : operatorIds )
   {
      localMatricesVector.push_back( face.getData( operatorId ) );
   }

   std::vector< real_t* > coeffs;
   for( auto coeffId : coeffIds )
   {
      coeffs.push_back( face.getData( coeffId )->getPointer( Level ) );
   }

   ValueType tmp;

   std::array< SD, 3 > triangleBlueSW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S};
   std::array< SD, 3 > triangleGrayS  = {SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE};
   std::array< SD, 3 > triangleBlueSE = {SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E};
   std::array< SD, 3 > triangleGrayNW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_NW};
   std::array< SD, 3 > triangleBlueN  = {SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N};
   std::array< SD, 3 > triangleGrayNE = {SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_E};

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         if( update == Replace )
         {
            tmp = ValueType( 0 );
         } else
         {
            tmp = dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];
         }

         for( uint_t coeffIdx = 0; coeffIdx < coeffIds.size(); ++coeffIdx )
         {
            tmp += assembleLocal< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix( Level ),
                                               src,
                                               coeffs[coeffIdx],
                                               triangleGrayS,
                                               {2, 0, 1} );
            tmp += assembleLocal< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix( Level ),
                                               src,
                                               coeffs[coeffIdx],
                                               triangleBlueSE,
                                               {1, 2, 0} );
            tmp += assembleLocal< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix( Level ),
                                               src,
                                               coeffs[coeffIdx],
                                               triangleBlueSW,
                                               {0, 1, 2} );
            tmp += assembleLocal< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix( Level ),
                                               src,
                                               coeffs[coeffIdx],
                                               triangleGrayNW,
                                               {1, 0, 2} );
            tmp += assembleLocal< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix( Level ),
                                               src,
                                               coeffs[coeffIdx],
                                               triangleBlueN,
                                               {2, 1, 0} );
            tmp += assembleLocal< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix( Level ),
                                               src,
                                               coeffs[coeffIdx],
                                               triangleGrayNE,
                                               {0, 2, 1} );
         }

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )] = tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void applyCoefficientDG( const uint_t&                                               Level,
                                Face&                                                       face,
                                const PrimitiveDataID< FaceP1LocalMatrixMemory, Face >&     operatorId,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                                const PrimitiveDataID< FunctionMemory< ValueType >, Face >& coeffId,
                                UpdateType                                                  update )
{
   typedef stencilDirection SD;

   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto localMatrices = face.getData( operatorId );
   auto src           = face.getData( srcId )->getPointer( Level );
   auto dst           = face.getData( dstId )->getPointer( Level );
   auto coeff         = face.getData( coeffId )->getPointer( Level );

   ValueType tmp;

   std::array< SD, 3 > triangleBlueSW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S};
   std::array< SD, 3 > triangleGraySE = {SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE};
   std::array< SD, 3 > triangleBlueSE = {SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E};
   std::array< SD, 3 > triangleGrayNW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_NW};
   std::array< SD, 3 > triangleBlueNW = {SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N};
   std::array< SD, 3 > triangleGrayNE = {SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_E};

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         if( update == Replace )
         {
            tmp = ValueType( 0 );
         } else
         {
            tmp = dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];
         }

         tmp +=
             coeff[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_GRAY_SE )] *
             assembleLocalDG< ValueType >( Level, i, j, localMatrices->getGrayMatrix( Level ), src, triangleGraySE, {2, 0, 1} );
         tmp +=
             coeff[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_BLUE_SE )] *
             assembleLocalDG< ValueType >( Level, i, j, localMatrices->getBlueMatrix( Level ), src, triangleBlueSE, {1, 2, 0} );
         tmp +=
             coeff[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_BLUE_SW )] *
             assembleLocalDG< ValueType >( Level, i, j, localMatrices->getBlueMatrix( Level ), src, triangleBlueSW, {0, 1, 2} );
         tmp +=
             coeff[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_GRAY_NW )] *
             assembleLocalDG< ValueType >( Level, i, j, localMatrices->getGrayMatrix( Level ), src, triangleGrayNW, {1, 0, 2} );
         tmp +=
             coeff[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_BLUE_NW )] *
             assembleLocalDG< ValueType >( Level, i, j, localMatrices->getBlueMatrix( Level ), src, triangleBlueNW, {2, 1, 0} );
         tmp +=
             coeff[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_GRAY_NE )] *
             assembleLocalDG< ValueType >( Level, i, j, localMatrices->getGrayMatrix( Level ), src, triangleGrayNE, {0, 2, 1} );

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )] = tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void smooth_gs( const uint_t&                                               Level,
                       Face&                                                       face,
                       const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto opr_data = face.getData( operatorId )->getPointer( Level );
   auto dst      = face.getData( dstId )->getPointer( Level );
   auto rhs      = face.getData( rhsId )->getPointer( Level );

   const auto invCenterWeight = 1.0 / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];

   ValueType tmp;

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         if( face.getNumNeighborCells() == 0 )
         {
            for( const auto direction : vertexdof::macroface::neighborsWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         } else if( face.getNumNeighborCells() == 1 )
         {
            for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         } else if( face.getNumNeighborCells() == 2 )
         {
            for( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
            {
               tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                      dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
            }
         }

         WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp * invCenterWeight;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void smooth_gs_coefficient( uint_t                                                                     Level,
                                   Face&                                                                      face,
                                   const std::vector< PrimitiveDataID< FaceP1LocalMatrixMemory, Face > >&     operatorIds,
                                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstId,
                                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                rhsId,
                                   const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& coeffIds )
{
   typedef stencilDirection SD;

   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto dst = face.getData( dstId )->getPointer( Level );
   auto rhs = face.getData( rhsId )->getPointer( Level );

   std::vector< FaceP1LocalMatrixMemory* > localMatricesVector;
   for( auto operatorId : operatorIds )
   {
      localMatricesVector.push_back( face.getData( operatorId ) );
   }

   std::vector< real_t* > coeffs;
   for( auto coeffId : coeffIds )
   {
      coeffs.push_back( face.getData( coeffId )->getPointer( Level ) );
   }

   std::array< SD, 3 > triangleBlueSW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_S};
   std::array< SD, 3 > triangleGrayS  = {SD::VERTEX_C, SD::VERTEX_S, SD::VERTEX_SE};
   std::array< SD, 3 > triangleBlueSE = {SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E};
   std::array< SD, 3 > triangleGrayNW = {SD::VERTEX_C, SD::VERTEX_W, SD::VERTEX_NW};
   std::array< SD, 3 > triangleBlueN  = {SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N};
   std::array< SD, 3 > triangleGrayNE = {SD::VERTEX_C, SD::VERTEX_N, SD::VERTEX_E};

   ValueType             tmp;
   std::vector< real_t > opr_data( 7 );

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         std::fill( opr_data.begin(), opr_data.end(), 0.0 );

         for( uint_t coeffIdx = 0; coeffIdx < coeffIds.size(); ++coeffIdx )
         {
            assembleLocalStencil< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix( Level ),
                                               opr_data.data(),
                                               coeffs[coeffIdx],
                                               triangleGrayS,
                                               {2, 0, 1} );
            assembleLocalStencil< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix( Level ),
                                               opr_data.data(),
                                               coeffs[coeffIdx],
                                               triangleBlueSE,
                                               {1, 2, 0} );
            assembleLocalStencil< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix( Level ),
                                               opr_data.data(),
                                               coeffs[coeffIdx],
                                               triangleBlueSW,
                                               {0, 1, 2} );
            assembleLocalStencil< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix( Level ),
                                               opr_data.data(),
                                               coeffs[coeffIdx],
                                               triangleGrayNW,
                                               {1, 0, 2} );
            assembleLocalStencil< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getBlueMatrix( Level ),
                                               opr_data.data(),
                                               coeffs[coeffIdx],
                                               triangleBlueN,
                                               {2, 1, 0} );
            assembleLocalStencil< ValueType >( Level,
                                               i,
                                               j,
                                               localMatricesVector[coeffIdx]->getGrayMatrix( Level ),
                                               opr_data.data(),
                                               coeffs[coeffIdx],
                                               triangleGrayNE,
                                               {0, 2, 1} );
         }

         tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         //for (auto neighbor : neighbors) {
         for( uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k )
         {
            tmp -= opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[k] )] *
                   dst[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[k] )];
         }

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] =
             tmp / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void smooth_sor( const uint_t&                                               Level,
                        Face&                                                       face,
                        const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                        ValueType                                                   relax )
{
  uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData( operatorId )->getPointer( Level );
  auto dst      = face.getData( dstId )->getPointer( Level );
  auto rhs      = face.getData( rhsId )->getPointer( Level );

  const auto invCenterWeight = 1.0 / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];

  ValueType tmp;

  for( uint_t j = 1; j < rowsize - 2; ++j )
  {
    for( uint_t i = 1; i < inner_rowsize - 2; ++i )
    {
      tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

      if( face.getNumNeighborCells() == 0 )
      {
        for( const auto direction : vertexdof::macroface::neighborsWithoutCenter )
        {
          tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                 dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
        }
      } else if( face.getNumNeighborCells() == 1 )
      {
        for( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
        {
          tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                 dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
        }
      } else if( face.getNumNeighborCells() == 2 )
      {
        for( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
        {
          tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                 dst[vertexdof::macroface::indexFromVertex( Level, i, j, direction )];
        }
      }

      WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

      dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] =
      ( 1.0 - relax ) * dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] +
      relax * tmp * invCenterWeight;
    }
    --inner_rowsize;
  }
}

template < typename ValueType >
inline void smooth_jac( const uint_t&                                               Level,
                        Face&                                                       face,
                        const PrimitiveDataID< StencilMemory< ValueType >, Face >&  operatorId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                        const PrimitiveDataID< FunctionMemory< ValueType >, Face >& tmpId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto opr_data = face.getData( operatorId )->getPointer( Level );
   auto dst      = face.getData( dstId )->getPointer( Level );
   auto rhs      = face.getData( rhsId )->getPointer( Level );
   auto tmpVar   = face.getData( tmpId )->getPointer( Level );

   ValueType tmp;

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         for( auto neighbor : vertexdof::macroface::neighborsWithoutCenter )
         {
            tmp -= opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
                   tmpVar[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
         }

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] =
             tmp / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
      }
      --inner_rowsize;
   }
}

/// Checks if a given index is a the boundary of the face
/// \param index The index which should be checked
/// \param length Size of the triangle in the first dimension
inline bool is_boundary( uint_t index, uint_t length )
{
   if( index < length )
      return true;
   while( index >= length )
   {
      index -= length;
      length--;
   }
   return ( index == 0 || index == ( length - 1 ) );
}

template < typename ValueType >
inline void
    enumerate( const uint_t& Level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId, ValueType& num )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   uint_t mr = 1 + rowsize;

   ValueType* dstPtr = face.getData( dstId )->getPointer( Level );

   for( uint_t i = 0; i < rowsize - 3; ++i )
   {
      for( uint_t j = 0; j < inner_rowsize - 3; ++j )
      {
         dstPtr[mr] = num;
         num++;

         mr += 1;
      }

      mr += 2;
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void integrateDG( const uint_t&                                               Level,
                         Face&                                                       face,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsId,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >& rhsP1Id,
                         const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   using namespace vertexdof::macroface;
   typedef stencilDirection SD;

   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto rhs   = face.getData( rhsId )->getPointer( Level );
   auto rhsP1 = face.getData( rhsP1Id )->getPointer( Level );
   auto dst   = face.getData( dstId )->getPointer( Level );

   real_t faceArea         = std::pow( 4.0, -walberla::real_c( Level ) ) * face.area;
   real_t weightedFaceArea = faceArea / 3.0;

   ValueType tmp;

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         tmp =
             rhs[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_BLUE_SW )] *
             ( 0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_W )] ) +
               0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                     rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_S )] ) );
         tmp +=
             rhs[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_GRAY_SE )] *
             ( 0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_S )] ) +
               0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                     rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_SE )] ) );
         tmp += rhs[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_BLUE_SE )] *
                ( 0.5 * 0.5 *
                      ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                        rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_SE )] ) +
                  0.5 * 0.5 *
                      ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                        rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_E )] ) );

         tmp +=
             rhs[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_GRAY_NW )] *
             ( 0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_W )] ) +
               0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                     rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_NW )] ) );
         tmp += rhs[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_BLUE_NW )] *
                ( 0.5 * 0.5 *
                      ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                        rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_NW )] ) +
                  0.5 * 0.5 *
                      ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                        rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_N )] ) );
         tmp +=
             rhs[facedof::macroface::indexFaceFromVertex( Level, i, j, SD::CELL_GRAY_NE )] *
             ( 0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_N )] ) +
               0.5 * 0.5 *
                   ( rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_C )] +
                     rhsP1[indexFromVertex( Level, i, j, SD::VERTEX_E )] ) );

         dst[indexFromVertex( Level, i, j, SD::VERTEX_C )] = weightedFaceArea * tmp;
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline ValueType getMaxValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   auto src      = face.getData( srcId )->getPointer( level );
   auto localMax = -std::numeric_limits< ValueType >::max();

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         localMax = std::max( localMax, src[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )]);
      }
      --inner_rowsize;
   }

   return localMax;
}

template < typename ValueType >
inline ValueType
    getMaxMagnitude( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   auto src      = face.getData( srcId )->getPointer( level );
   auto localMax = ValueType( 0.0 );

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         localMax = std::max( localMax,
                              std::abs( src[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] ) );
      }
      --inner_rowsize;
   }

   return localMax;
}

template < typename ValueType >
inline ValueType getMinValue( const uint_t& level, Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   uint_t inner_rowsize = rowsize;

   auto src      = face.getData( srcId )->getPointer( level );
   auto localMin = std::numeric_limits< ValueType >::max();

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         localMin = std::min( localMin, src[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] );
      }
      --inner_rowsize;
   }

   return localMin;
}

#ifdef HHG_BUILD_WITH_PETSC


inline void saveOperator( const uint_t&                                              Level,
                          Face&                                                      face,
                          const PrimitiveDataID< StencilMemory< real_t >, Face >&    operatorId,
                          const PrimitiveDataID< FunctionMemory< PetscInt >, Face >& srcId,
                          const PrimitiveDataID< FunctionMemory< PetscInt >, Face >& dstId,
                          Mat&                                                       mat )
{
  uint_t rowsize = levelinfo::num_microvertices_per_edge( Level );
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData( operatorId )->getPointer( Level );
  auto src = face.getData( srcId )->getPointer( Level );
  auto dst = face.getData( dstId )->getPointer( Level );

  for ( uint_t i = 1; i < rowsize - 2; ++i )
  {
    for ( uint_t j = 1; j < inner_rowsize - 2; ++j )
    {

      PetscInt srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
      PetscInt dstInt = dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, VERTEX_C)], opr_data[VERTEX_C]);
      MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )], ADD_VALUES );

      if ( face.getNumNeighborCells() == 0 )
      {
        for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter )
        {
          srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
          //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
          MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[vertexdof::stencilIndexFromVertex( neighbor )], ADD_VALUES );
        }
      } else if ( face.getNumNeighborCells() == 1 )
      {
        for ( const auto neighbor : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
        {
          srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
          //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
          MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[vertexdof::stencilIndexFromVertex( neighbor )], ADD_VALUES );
        }
      } else if ( face.getNumNeighborCells() == 2 )
      {
        for ( const auto neighbor : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
        {
          srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor )];
          //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
          MatSetValues( mat, 1, &dstInt, 1, &srcInt, &opr_data[vertexdof::stencilIndexFromVertex( neighbor )], ADD_VALUES );
        }
      }
    }
    --inner_rowsize;
  }
}

template < typename ValueType >
inline void createVectorFromFunction( const uint_t&                                               Level,
                                      Face&                                                       face,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                                      const PrimitiveDataID< FunctionMemory< PetscInt >, Face >&  numeratorId,
                                      Vec&                                                        vec )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto src       = face.getData( srcId )->getPointer( Level );
   auto numerator = face.getData( numeratorId )->getPointer( Level );

   for( uint_t i = 1; i < rowsize - 2; ++i )
   {
      for( uint_t j = 1; j < inner_rowsize - 2; ++j )
      {
         PetscInt numeratorInt = numerator[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         VecSetValues( vec,
                       1,
                       &numeratorInt,
                       &src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )],
                       INSERT_VALUES );
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void createFunctionFromVector( const uint_t&                                               Level,
                                      Face&                                                       face,
                                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                                      const PrimitiveDataID< FunctionMemory< PetscInt >, Face >&  numeratorId,
                                      Vec&                                                        vec )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto src       = face.getData( srcId )->getPointer( Level );
   auto numerator = face.getData( numeratorId )->getPointer( Level );

   for( uint_t i = 1; i < rowsize - 2; ++i )
   {
      for( uint_t j = 1; j < inner_rowsize - 2; ++j )
      {
         PetscInt numeratorInt = numerator[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
         VecGetValues(
             vec, 1, &numeratorInt, &src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] );
      }
      --inner_rowsize;
   }
}

inline void applyDirichletBC(const uint_t & level, Face & face,std::vector<PetscInt> &mat,
                             const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &numeratorId){

  for( const auto & it : vertexdof::macroface::Iterator( level, 1 ))
  {
    mat.push_back(face.getData(numeratorId)->getPointer( level )[vertexdof::macroface::index( level, it.x(), it.y() )]);
  }

}

#endif

template < typename ValueType >
inline void
    printFunctionMemory( const uint_t Level, const Face& face, const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId )
{
   ValueType* faceMemory = face.getData( dstId )->getPointer( Level );
   using namespace std;
   cout << setfill( '=' ) << setw( 100 ) << "" << endl;
   cout << face << std::left << setprecision( 1 ) << fixed << setfill( ' ' ) << endl << "Vertex DoFs: ";
   for( const auto& it : vertexdof::macroface::Iterator( Level, 0 ) )
   {
      if( it.col() == 0 )
         std::cout << std::endl;
      cout << setw( 5 )
           << faceMemory[hhg::vertexdof::macroface::indexFromVertex( Level, it.col(), it.row(), stencilDirection::VERTEX_C )]
           << "|";
   }
   cout << endl << setfill( '=' ) << setw( 100 ) << "" << endl << setfill( ' ' );
}

} // namespace macroface
} // namespace vertexdof
} // namespace hhg
