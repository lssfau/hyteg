#pragma once

#include "core/debug/all.h"
#include "core/math/KahanSummation.h"

#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/P1Elements.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFaceIndex.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"
#include "tinyhhg_core/indexing/Common.hpp"

namespace hhg {
namespace vertexdof {
namespace macroface {

using walberla::uint_t;
using walberla::real_c;
using indexing::Index;

inline Point3D coordinateFromIndex( const uint_t & Level, const Face & face, const Index & index )
{
  const real_t  stepFrequency = 1.0 / levelinfo::num_microedges_per_edge( Level );
  const Point3D xStep         = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) * stepFrequency;
  const Point3D yStep         = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) * stepFrequency;
  return face.getCoordinates()[0] + xStep * real_c( index.x() ) + yStep * real_c( index.y() );
}

template< typename ValueType >
inline ValueType assembleLocal(const uint_t & Level,
                               uint_t i, uint_t j, const Matrix3r& localMatrix,
                               double* src,
                               double* coeff,
                               const std::array<stencilDirection,3>& vertices,
                               const std::array<uint_t,3>& idx)
{
  ValueType meanCoeff = 1.0/3.0 * (coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )]
                                 + coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )]
                                 + coeff[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )]);

  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )]
         + localMatrix(idx[0],idx[1]) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )]
         + localMatrix(idx[0],idx[2]) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )];
  return meanCoeff * tmp;
}

template< typename ValueType >
inline ValueType assembleLocalDG( const uint_t & Level, uint_t i, uint_t j, const Matrix3r& localMatrix,
                               double* src,
                               const std::array<stencilDirection,3>& vertices,
                               const std::array<uint_t,3>& idx)
{
  ValueType tmp;
  tmp  = localMatrix(idx[0],idx[0]) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[0] )]
      + localMatrix(idx[0],idx[1]) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[1] )]
      + localMatrix(idx[0],idx[2]) * src[vertexdof::macroface::indexFromVertex( Level, i, j, vertices[2] )];
  return tmp;
}

template< typename ValueType >
inline void interpolate(const uint_t & Level,
                            Face &face,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face>& faceMemoryId,
                            const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                            std::function<ValueType(const hhg::Point3D &, const std::vector<ValueType>&)> &expr)
{
  ValueType * faceData = face.getData( faceMemoryId )->getPointer( Level );

  std::vector< ValueType * > srcPtr;
  for( const auto & src : srcIds )
  {
    srcPtr.push_back( face.getData(src)->getPointer( Level ) );
  }

  std::vector<ValueType> srcVector( srcIds.size() );

  Point3D xBlend;

  for ( const auto & it : vertexdof::macroface::Iterator( Level, 1 ) )
  {
    const Point3D coordinate = coordinateFromIndex( Level, face, it );
    const uint_t  idx        = vertexdof::macroface::indexFromVertex( Level, it.x(), it.y(),
                                                                               stencilDirection::VERTEX_C );

    for ( uint_t k = 0; k < srcPtr.size(); ++k )
    {
      srcVector[ k ] = srcPtr[ k ][ idx ];
    }
    face.getGeometryMap()->evalF(coordinate, xBlend);
    faceData[ idx ] = expr( xBlend, srcVector );
  }
}


template< typename ValueType >
inline void assign( const uint_t & Level, Face &face,
                   const std::vector<ValueType>& scalars,
                   const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                   const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* dst = face.getData(dstId)->getPointer( Level );
  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }
  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      ValueType tmp = scalars[0]*srcPtr[0][vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                                           stencilDirection::VERTEX_C )];

      for (uint_t k = 1; k < srcIds.size(); ++k) {
        tmp += scalars[k]*srcPtr[k][vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                                    stencilDirection::VERTEX_C )];
      }
      dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
    }
    --inner_rowsize;
  }
}


template< typename ValueType >
inline void add(const uint_t & Level, Face &face,
                const std::vector<ValueType>& scalars,
                const std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>> &srcIds,
                const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* dstPtr = face.getData(dstId)->getPointer( Level );
  std::vector<ValueType*> srcPtr;
  for(auto src : srcIds){
    srcPtr.push_back(face.getData(src)->getPointer( Level ));
  }

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      ValueType tmp = 0.0;

      for (uint_t k = 0; k < srcIds.size(); ++k) {
        tmp += scalars[k] * srcPtr[k][vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                                      stencilDirection::VERTEX_C )];
      }

      dstPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
    }

    --inner_rowsize;
  }
}


template< typename ValueType >
inline real_t dot(const uint_t & Level, Face &face,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Face>& lhsId,
                  const PrimitiveDataID<FunctionMemory< ValueType >, Face>& rhsId) {

  walberla::math::KahanAccumulator< ValueType > scalarProduct;
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* lhsPtr = face.getData(lhsId)->getPointer( Level );
  ValueType* rhsPtr = face.getData(rhsId)->getPointer( Level );

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      scalarProduct += lhsPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )]
          * rhsPtr[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];
    }
    --inner_rowsize;
  }

  return scalarProduct.get();
}


template< typename ValueType >
inline void apply(const uint_t & Level, Face &face, const PrimitiveDataID< StencilMemory< ValueType >, Face>& operatorId,
                       const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                       const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId, UpdateType update) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  ValueType* opr_data = face.getData(operatorId)->getPointer( Level );
  ValueType* src = face.getData(srcId)->getPointer( Level );
  ValueType* dst = face.getData(dstId)->getPointer( Level );

  ValueType tmp;

  if( update == Replace ) {
    for (uint_t j = 1; j < rowsize - 2; ++j) {
      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
        tmp = opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

        //strangely the intel compiler cant handle this if it is a loop
        static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
          tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[0] )]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                    vertexdof::macroface::neighborsWithoutCenter[0] )];
          tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[1] )]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                    vertexdof::macroface::neighborsWithoutCenter[1] )];
          tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[2] )]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                    vertexdof::macroface::neighborsWithoutCenter[2] )];
          tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[3] )]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                    vertexdof::macroface::neighborsWithoutCenter[3] )];
          tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[4] )]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                    vertexdof::macroface::neighborsWithoutCenter[4] )];
          tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[5] )]
              * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                    vertexdof::macroface::neighborsWithoutCenter[5] )];

        if ( face.getNumNeighborCells() == 1 )
        {
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BC )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BC )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BS )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BS )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BSW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BSW )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BW )];
        }
        else if ( face.getNumNeighborCells() == 2 )
        {
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BC )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BC )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BS )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BS )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BSW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BSW )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BW )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FC )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FC )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FN )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FN )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FNE )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FNE )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FE )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FE )];
        }

        WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

        dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp;
      }
      --inner_rowsize;
    }
  } else {
    for (uint_t j = 1; j < rowsize - 2; ++j) {
      for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
        tmp = opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] * src[vertexdof::macroface::indexFromVertex(
        Level, i, j, stencilDirection::VERTEX_C )];

        //strangely the intel compiler cant handle this if it is a loop
        static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[0])]*src[vertexdof::macroface::indexFromVertex(
        Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0] )];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[1])]*src[vertexdof::macroface::indexFromVertex(
        Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1] )];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[2])]*src[vertexdof::macroface::indexFromVertex(
        Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2] )];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[3])]*src[vertexdof::macroface::indexFromVertex(
        Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3] )];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[4])]*src[vertexdof::macroface::indexFromVertex(
        Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4] )];
        tmp += opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[5])]*src[vertexdof::macroface::indexFromVertex(
        Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5] )];

        if ( face.getNumNeighborCells() == 1 )
        {
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BC )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BC )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BS )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BS )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BSW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BSW )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BW )];
        }
        else if ( face.getNumNeighborCells() == 2 )
        {
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BC )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BC )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BS )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BS )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BSW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BSW )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_BW )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_BW )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FC )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FC )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FN )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FN )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FNE )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FNE )];
            tmp += opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_FE )]
                * src[vertexdof::macroface::indexFromVertex( Level, i, j,
                                                                      stencilDirection::VERTEX_FE )];
        }

        WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

        dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] += tmp;
      }
      --inner_rowsize;
    }
  }
}


template< typename ValueType >
inline void applyCoefficient(const uint_t & Level, Face &face, const PrimitiveDataID<FaceP1LocalMatrixMemory, Face>& operatorId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &coeffId,
                                 UpdateType update) {
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto localMatrices = face.getData(operatorId);
  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);
  auto coeff = face.getData(coeffId)->getPointer(Level);

  ValueType tmp;

  std::array<SD,3> triangleBlueSW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_S  };
  std::array<SD,3> triangleGrayS  = { SD::VERTEX_C, SD::VERTEX_S,  SD::VERTEX_SE };
  std::array<SD,3> triangleBlueSE = { SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E  };
  std::array<SD,3> triangleGrayNW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_NW };
  std::array<SD,3> triangleBlueN  = { SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N  };
  std::array<SD,3> triangleGrayNE = { SD::VERTEX_C, SD::VERTEX_N,  SD::VERTEX_E  };

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {


      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];
      }

      tmp += assembleLocal< ValueType >( Level, i, j, localMatrices->getGrayMatrix(Level), src, coeff, triangleGrayS, {2,0,1});
      tmp += assembleLocal< ValueType >( Level, i, j, localMatrices->getBlueMatrix(Level), src, coeff, triangleBlueSE, {1,2,0});
      tmp += assembleLocal< ValueType >( Level, i, j, localMatrices->getBlueMatrix(Level), src, coeff, triangleBlueSW, {0,1,2});
      tmp += assembleLocal< ValueType >( Level, i, j, localMatrices->getGrayMatrix(Level), src, coeff, triangleGrayNW, {1,0,2});
      tmp += assembleLocal< ValueType >( Level, i, j, localMatrices->getBlueMatrix(Level), src, coeff, triangleBlueN, {2,1,0});
      tmp += assembleLocal< ValueType >( Level, i, j, localMatrices->getGrayMatrix(Level), src, coeff, triangleGrayNE, {0,2,1});

      dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )] = tmp;
    }
    --inner_rowsize;
  }
}


template< typename ValueType >
inline void applyCoefficientDG(const uint_t & Level, Face &face, const PrimitiveDataID<FaceP1LocalMatrixMemory, Face>& operatorId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                                 const PrimitiveDataID<FunctionMemory< ValueType >, Face> &coeffId,
                                 UpdateType update) {
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto localMatrices = face.getData(operatorId);
  auto src = face.getData(srcId)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);
  auto coeff = face.getData(coeffId)->getPointer(Level);

  ValueType tmp;

  std::array<SD,3> triangleBlueSW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_S  };
  std::array<SD,3> triangleGraySE = { SD::VERTEX_C, SD::VERTEX_S,  SD::VERTEX_SE };
  std::array<SD,3> triangleBlueSE = { SD::VERTEX_C, SD::VERTEX_SE, SD::VERTEX_E  };
  std::array<SD,3> triangleGrayNW = { SD::VERTEX_C, SD::VERTEX_W,  SD::VERTEX_NW };
  std::array<SD,3> triangleBlueNW = { SD::VERTEX_C, SD::VERTEX_NW, SD::VERTEX_N  };
  std::array<SD,3> triangleGrayNE = { SD::VERTEX_C, SD::VERTEX_N,  SD::VERTEX_E  };

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {


      if (update == Replace) {
        tmp = ValueType(0);
      }
      else {
        tmp = dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];
      }

      tmp += coeff[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_GRAY_SE)] * assembleLocalDG< ValueType >( Level, i, j, localMatrices->getGrayMatrix(Level), src, triangleGraySE, {2,0,1});
      tmp += coeff[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_BLUE_SE)] * assembleLocalDG< ValueType >( Level, i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueSE, {1,2,0});
      tmp += coeff[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_BLUE_SW)] * assembleLocalDG< ValueType >( Level, i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueSW, {0,1,2});
      tmp += coeff[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_GRAY_NW)] * assembleLocalDG< ValueType >( Level, i, j, localMatrices->getGrayMatrix(Level), src, triangleGrayNW, {1,0,2});
      tmp += coeff[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_BLUE_NW)] * assembleLocalDG< ValueType >( Level, i, j, localMatrices->getBlueMatrix(Level), src, triangleBlueNW, {2,1,0});
      tmp += coeff[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_GRAY_NE)] * assembleLocalDG< ValueType >( Level, i, j, localMatrices->getGrayMatrix(Level), src, triangleGrayNE, {0,2,1});

      dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )] = tmp;
    }
    --inner_rowsize;
  }

}


template< typename ValueType >
inline void smooth_gs(const uint_t & Level, Face &face, const PrimitiveDataID<StencilMemory< ValueType >, Face>& operatorId,
                           const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                           const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]*dst[vertexdof::macroface::indexFromVertex(
         Level, i, j, vertexdof::macroface::neighborsWithoutCenter[k] )];
      }

      dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}


template< typename ValueType >
inline void smooth_sor(const uint_t & Level, Face &face, const PrimitiveDataID<StencilMemory< ValueType >, Face>& operatorId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                            ValueType relax) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

      //for (auto neighbor : neighbors) {
      for(uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k){
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])]*dst[vertexdof::macroface::indexFromVertex(
         Level, i, j, vertexdof::macroface::neighborsWithoutCenter[k] )];
      }

      dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = ( 1.0 - relax) * dst[vertexdof::macroface::indexFromVertex(
       Level, i, j, stencilDirection::VERTEX_C )] + relax * tmp / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}


template< typename ValueType >
inline void smooth_jac(const uint_t & Level, Face &face, const PrimitiveDataID<StencilMemory< ValueType >, Face>& operatorId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &tmpId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );
  auto rhs = face.getData(rhsId)->getPointer( Level );
  auto tmpVar = face.getData(tmpId)->getPointer( Level );

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

      for (auto neighbor : vertexdof::macroface::neighborsWithoutCenter ) {
        tmp -= opr_data[vertexdof::stencilIndexFromVertex(neighbor)]*tmpVar[vertexdof::macroface::indexFromVertex(
        Level, i, j, neighbor )];
      }

      dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] = tmp / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C)];
    }
    --inner_rowsize;
  }
}


template< typename ValueType >
inline void prolongate(const uint_t & Level, Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face>& memoryId) {

  typedef stencilDirection SD;
  using namespace vertexdof::macroface;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;

  auto v_f = face.getData(memoryId)->getPointer( Level + 1 );
  auto v_c = face.getData(memoryId)->getPointer( Level );

  uint_t j;

  for (uint_t i = 1; i < N_c - 1; ++i) {
    for (j = 1; j < N_c_i - 2; ++j) {
      v_f[indexFromVertex( Level + 1, 2 * i, 2 * j, SD::VERTEX_C )] = v_c[indexFromVertex( Level, i, j, SD::VERTEX_C )];
      v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] =
          0.5*( v_c[indexFromVertex( Level, i - 1, j, SD::VERTEX_C )] + v_c[indexFromVertex( Level, i, j - 1, SD::VERTEX_C )]);
      v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex( Level, i, j, SD::VERTEX_C )]
                                                                                  + v_c[indexFromVertex( Level, i - 1, j, SD::VERTEX_C )]);
      v_f[indexFromVertex( Level + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex(
      Level,
      i, j,
      SD::VERTEX_C )] + v_c[indexFromVertex(
      Level, i, j - 1, SD::VERTEX_C )]);
    }

    v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex(
    Level, i - 1, j, SD::VERTEX_C )] + v_c[indexFromVertex(
    Level, i, j - 1, SD::VERTEX_C )]);
    v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex(
    Level, i,
    j, SD::VERTEX_C )] + v_c[indexFromVertex(
    Level, i - 1, j, SD::VERTEX_C )]);
    v_f[indexFromVertex( Level + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex(
    Level, i,
    j, SD::VERTEX_C )] + v_c[indexFromVertex(
    Level, i, j - 1, SD::VERTEX_C )]);

    --N_c_i;
  }
}


template< typename ValueType >
inline void prolongateQuadratic(const uint_t & Level, Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face>& memoryId) {

  typedef stencilDirection SD;
  using namespace vertexdof::macroface;

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level);
  uint_t N_c_i = N_c;
  auto v_f = face.getData(memoryId)->getPointer( Level + 1 );
  auto v_c = face.getData(memoryId)->getPointer( Level );

  uint_t i, j;
  ValueType linearx, lineary, linearxy, offx, offy, offxy;
  i = 0;
  for ( j = 2; j <= N_c - 1; j += 2 ) {
// upper triangle inner points
//calculate offsets
    linearx = 0.5 * ( v_c[indexFromVertex( Level, i, j - 2, SD::VERTEX_C )] + v_c[indexFromVertex( Level, i, j, SD::VERTEX_C )] );
    lineary = 0.5 * ( v_c[indexFromVertex( Level, i + 2, j - 2, SD::VERTEX_C )] + v_c[indexFromVertex( Level,
                                                                                                       i, j - 2, SD::VERTEX_C )] );
    linearxy = 0.5 * ( v_c[indexFromVertex( Level, i + 2, j - 2, SD::VERTEX_C )] + v_c[indexFromVertex( Level,
                                                                                                        i, j, SD::VERTEX_C )] );

    offx = v_c[indexFromVertex( Level, i, j - 1, SD::VERTEX_C )] - linearx;
    offy = v_c[indexFromVertex( Level, i + 1, j - 2, SD::VERTEX_C )] - lineary;
    offxy = v_c[indexFromVertex( Level, i + 1, j - 1, SD::VERTEX_C )] - linearxy;

// left bottom corner
    v_f[indexFromVertex( Level + 1, 2 * i + 1, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearx + lineary ) + 0.5 * offx + 0.5 * offy + 0.25 * offxy;
// right bottom corner
    v_f[indexFromVertex( Level + 1, 2 * i + 1, 2 * j - 2, SD::VERTEX_C )] = 0.5 * ( linearx + linearxy ) + 0.5 * offx + 0.25 * offy + 0.5 * offxy;
// top corner
    v_f[indexFromVertex( Level + 1, 2 * i + 2, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearxy + lineary ) + 0.25 * offx + 0.5 * offy + 0.5 * offxy;
  }

  N_c_i -= 1;

  for ( j = 2; j < N_c - 1; j += 2 ) {
    for ( i = 2; i < N_c_i - 1; i += 2 ) {

// upper triangle inner points
//calculate offsets
      linearx = 0.5 * ( v_c[indexFromVertex( Level, i, j - 2, SD::VERTEX_C )] + v_c[indexFromVertex( Level,
                                                                                                     i, j, SD::VERTEX_C )] );
      lineary = 0.5 * ( v_c[indexFromVertex( Level, i + 2, j - 2, SD::VERTEX_C )] + v_c[indexFromVertex( Level,
                                                                                                         i, j - 2, SD::VERTEX_C )] );
      linearxy = 0.5 * ( v_c[indexFromVertex( Level, i + 2, j - 2, SD::VERTEX_C )] + v_c[indexFromVertex( Level,
                                                                                                          i, j, SD::VERTEX_C )] );

      offx = v_c[indexFromVertex( Level, i, j - 1, SD::VERTEX_C )] - linearx;
      offy = v_c[indexFromVertex( Level, i + 1, j - 2, SD::VERTEX_C )] - lineary;
      offxy = v_c[indexFromVertex( Level, i + 1, j - 1, SD::VERTEX_C )] - linearxy;
// left bottom corner
      v_f[indexFromVertex( Level + 1, 2 * i + 1, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearx + lineary ) + 0.5 * offx + 0.5 * offy + 0.25 * offxy;
// right bottom corner
      v_f[indexFromVertex( Level + 1, 2 * i + 1, 2 * j - 2, SD::VERTEX_C )] = 0.5 * ( linearx + linearxy ) + 0.5 * offx + 0.25 * offy + 0.5 * offxy;
// top corner
      v_f[indexFromVertex( Level + 1, 2 * i + 2, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearxy + lineary ) + 0.25 * offx + 0.5 * offy + 0.5 * offxy;

// lower triangle all points
//calculate offsets
      lineary = 0.5 * ( v_c[indexFromVertex( Level, i - 2, j, SD::VERTEX_C )] + v_c[indexFromVertex( Level,
                                                                                                     i, j, SD::VERTEX_C )] );
      linearxy = 0.5 * ( v_c[indexFromVertex( Level, i - 2, j, SD::VERTEX_C )] + v_c[indexFromVertex( Level,
                                                                                                      i, j - 2, SD::VERTEX_C )] );

      offy = v_c[indexFromVertex( Level, i - 1, j, SD::VERTEX_C )] - lineary;
      offxy = v_c[indexFromVertex( Level, i - 1, j - 1, SD::VERTEX_C )] - linearxy;
// first inner points
// left bottom corner
      v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( linearx + lineary ) + 0.5 * offx + 0.5 * offy + 0.25 * offxy;
// right bottom corner
      v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j - 2, SD::VERTEX_C )] = 0.5 * ( linearx + linearxy ) + 0.5 * offx + 0.25 * offy + 0.5 * offxy;
// top corner
      v_f[indexFromVertex( Level + 1, 2 * i - 2, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( linearxy + lineary ) + 0.25 * offx + 0.5 * offy + 0.5 * offxy;

// boundary points
// x-direction
      v_f[indexFromVertex( Level + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( linearx + v_c[indexFromVertex( Level,
                                                                                                                 i, j, SD::VERTEX_C )] ) + 0.75 * offx;
      v_f[indexFromVertex( Level + 1, 2 * i, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearx + v_c[indexFromVertex( Level,
                                                                                                                 i, j - 2, SD::VERTEX_C )] ) + 0.75 * offx;
//y-direction
      v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex( Level,

                                                                                                       i, j,
                                                                                                       SD::VERTEX_C )] + lineary ) + 0.75 * offy;
      v_f[indexFromVertex( Level + 1, 2 * i - 3, 2 * j, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex( Level,

                                                                                                       i - 2, j,
                                                                                                       SD::VERTEX_C )] + lineary ) + 0.75 * offy;
//xy-direction
      v_f[indexFromVertex( Level + 1, 2 * i - 1, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex( Level,
                                                                                                           i, j - 2, SD::VERTEX_C )] + linearxy ) + 0.75 * offxy;
      v_f[indexFromVertex( Level + 1, 2 * i - 3, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( v_c[indexFromVertex( Level,
                                                                                                           i - 2, j, SD::VERTEX_C )] + linearxy ) + 0.75 * offxy;
// coarse points
      v_f[indexFromVertex( Level + 1, 2 * i, 2 * j, SD::VERTEX_C )] = v_c[indexFromVertex( Level,
                                                                                           i, j,
                                                                                           SD::VERTEX_C )];
      v_f[indexFromVertex( Level + 1, 2 * i, 2 * j - 2, SD::VERTEX_C )] = v_c[indexFromVertex( Level,
                                                                                               i,
                                                                                               j - 1, SD::VERTEX_C )];
      v_f[indexFromVertex( Level + 1, 2 * i - 2, 2 * j, SD::VERTEX_C )] = v_c[indexFromVertex( Level,

                                                                                               i - 1, j,
                                                                                               SD::VERTEX_C )];
      v_f[indexFromVertex( Level + 1, 2 * i - 2, 2 * j - 2, SD::VERTEX_C )] = v_c[indexFromVertex( Level,

                                                                                                   i - 1, j - 1,
                                                                                                   SD::VERTEX_C )];
    }
    N_c_i -= 2;

  }
}


template< typename ValueType >
inline void restrict(const uint_t & Level, Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &memoryId) {

  uint_t N_c = levelinfo::num_microvertices_per_edge(Level - 1);
  uint_t N_c_i = N_c;

  auto v_f = face.getData(memoryId)->getPointer( Level );
  auto v_c = face.getData(memoryId)->getPointer( Level - 1 );

  ValueType tmp;

  for (uint_t j = 1; j < N_c - 2; ++j) {
    for (uint_t i = 1; i < N_c_i - 2; ++i) {

      tmp = v_f[vertexdof::macroface::indexFromVertex( Level, 2 * i, 2 * j,
                                                                stencilDirection::VERTEX_C )];

      for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter ) {
        tmp += 0.5*v_f[vertexdof::macroface::indexFromVertex( Level, 2 * i, 2 * j, neighbor )];
      }

      v_c[vertexdof::macroface::indexFromVertex( Level - 1, i, j, stencilDirection::VERTEX_C )] = tmp;
    }

    --N_c_i;
  }
}


/// Checks if a given index is a the boundary of the face
/// \param index The index which should be checked
/// \param length Size of the triangle in the first dimension
inline bool is_boundary(uint_t index, uint_t length) {
  if (index < length) return true;
  while (index >= length) {
    index -= length;
    length--;
  }
  return (index==0 || index==(length - 1));
}

template< typename ValueType >
inline void enumerate(const uint_t & Level, Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId, uint_t& num) {
  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  uint_t mr = 1 + rowsize;

  ValueType* dstPtr = face.getData(dstId)->getPointer( Level );

  for (uint_t i = 0; i < rowsize - 3; ++i) {
    for (uint_t j = 0; j < inner_rowsize - 3; ++j) {

      dstPtr[mr] = static_cast<ValueType>(num);
      num++;

      mr += 1;
    }

    mr += 2;
    --inner_rowsize;
  }
}


template< typename ValueType >
inline void integrateDG(const uint_t & Level, Face &face,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsId,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &rhsP1Id,
                            const PrimitiveDataID<FunctionMemory< ValueType >, Face> &dstId) {
  using namespace vertexdof::macroface;
  typedef stencilDirection SD;

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto rhs = face.getData(rhsId)->getPointer(Level);
  auto rhsP1 = face.getData(rhsP1Id)->getPointer(Level);
  auto dst = face.getData(dstId)->getPointer(Level);

  real_t faceArea = std::pow(4.0, -walberla::real_c(Level)) * face.area;
  real_t weightedFaceArea = faceArea / 3.0;

  ValueType tmp;

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {

      tmp  = rhs[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_BLUE_SW)] * (0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_W )]) + 0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_S )]));
      tmp += rhs[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_GRAY_SE)] * (0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_S )]) + 0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_SE )]));
      tmp += rhs[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_BLUE_SE)] * (0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_SE )]) + 0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_E )]));

      tmp += rhs[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_GRAY_NW)] * (0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_W )]) + 0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_NW )]));
      tmp += rhs[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_BLUE_NW)] * (0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_NW )]) + 0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_N )]));
      tmp += rhs[DGFace::indexDGFaceFromVertex( Level, i, j, SD::CELL_GRAY_NE)] * (0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_N )]) + 0.5 * 0.5 * ( rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_C )] + rhsP1[indexFromVertex( Level,
       i, j, SD::VERTEX_E )]));

      dst[indexFromVertex( Level,   i, j, SD::VERTEX_C )] = weightedFaceArea * tmp;
    }
    --inner_rowsize;
  }

}


template< typename ValueType >
inline real_t getMaxValue(const uint_t & Level, Face &face, const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  real_t localMax = std::numeric_limits<real_t>::min();

  for (uint_t j = 1; j < rowsize - 2; ++j) {
    for (uint_t i = 1; i < inner_rowsize - 2; ++i) {
      localMax = std::max(localMax, src[vertexdof::macroface::indexFromVertex( Level,   i, j,
                                                                                        stencilDirection::VERTEX_C )]);
    }
    --inner_rowsize;
  }

  return localMax;
}


#ifdef HHG_BUILD_WITH_PETSC

inline void saveOperator(const uint_t & Level, Face &face, const PrimitiveDataID<StencilMemory< real_t >, Face>& operatorId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &srcId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &dstId, Mat& mat) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto opr_data = face.getData(operatorId)->getPointer( Level );
  auto src = face.getData(srcId)->getPointer( Level );
  auto dst = face.getData(dstId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C)];
      PetscInt dstInt = dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C)];
      //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, VERTEX_C)], opr_data[VERTEX_C]);
      MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)] ,INSERT_VALUES);

      for ( const auto & neighbor : vertexdof::macroface::neighborsWithoutCenter ) {
        srcInt = src[vertexdof::macroface::indexFromVertex( Level, i, j, neighbor)];
        //out << fmt::format("{}\t{}\t{}\n", dst[index<Level>(i, j, VERTEX_C)], src[index<Level>(i, j, neighbor)], opr_data[neighbor]);
        MatSetValues(mat,1,&dstInt,1,&srcInt,&opr_data[vertexdof::stencilIndexFromVertex(neighbor)] ,INSERT_VALUES);
      }
    }
    --inner_rowsize;
  }
}



template< typename ValueType >
inline void createVectorFromFunction(const uint_t & Level, Face &face,
                              const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                              const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &numeratorId,
                              Vec& vec) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt numeratorInt = numerator[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C)];
      VecSetValues(vec,1,&numeratorInt,&src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C)],INSERT_VALUES);
    }
    --inner_rowsize;
  }
}


template< typename ValueType >
inline void createFunctionFromVector(const uint_t & Level, Face &face,
                                         const PrimitiveDataID<FunctionMemory< ValueType >, Face> &srcId,
                                         const PrimitiveDataID<FunctionMemory< PetscInt >, Face> &numeratorId,
                                         Vec& vec) {

  uint_t rowsize = levelinfo::num_microvertices_per_edge(Level);
  uint_t inner_rowsize = rowsize;

  auto src = face.getData(srcId)->getPointer( Level );
  auto numerator = face.getData(numeratorId)->getPointer( Level );


  for (uint_t i = 1; i < rowsize - 2; ++i) {
    for (uint_t j = 1; j < inner_rowsize - 2; ++j) {
      PetscInt numeratorInt = numerator[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C)];
      VecGetValues(vec,1,&numeratorInt,&src[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C)]);
    }
    --inner_rowsize;
  }
}

#endif

template< typename ValueType >
inline void printFunctionMemory(const uint_t Level,const Face &face, const PrimitiveDataID <FunctionMemory<ValueType>, Face> &dstId) {
  ValueType* faceMemory = face.getData(dstId)->getPointer( Level );
  using namespace std;
  cout << setfill('=') << setw(100) << "" << endl;
  cout << face << std::left << setprecision(1) << fixed << setfill(' ')  << endl << "Vertex DoFs: ";
  for ( const auto & it : vertexdof::macroface::Iterator( Level, 0 ) ){
    if(it.col() == 0) std::cout << std::endl;
    cout << setw(5) << faceMemory[hhg::vertexdof::macroface::indexFromVertex( Level,   it.col(),
                                                                                       it.row(),
                                                                                       stencilDirection::VERTEX_C )] << "|";
  }
  cout << endl << setfill('=') << setw(100) << "" << endl << setfill(' ');

}

}// namespace macroface
}// namespace vertexdof
}// namespace hhg

