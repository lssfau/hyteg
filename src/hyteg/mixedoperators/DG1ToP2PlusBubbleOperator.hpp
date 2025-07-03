/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/forms/form_hyteg_generated/dg1_to_p2_plus_bubble/dg1_to_p2_plus_bubble_divt_affine_q3.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2PlusBubbleFunction.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

namespace hyteg {

using namespace dg;
using facedof::FaceType;
using indexing::Index;
using volumedofspace::indexing::VolumeDoFMemoryLayout;
using walberla::int_c;
using walberla::real_t;

template < typename Form >
class DG1ToP2PlusBubbleOperator : public Operator< DG1Function< real_t >, P2PlusBubbleFunction< real_t > >
{
 public:
   DG1ToP2PlusBubbleOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                              uint_t                                     minLevel,
                              uint_t                                     maxLevel,
                              std::shared_ptr< Form >                    form = std::make_shared< Form >() )
   : Operator< DG1Function< real_t >, P2PlusBubbleFunction< real_t > >( storage, minLevel, maxLevel )
   , form_( form )
   {}

   typedef Form FormType;

   void apply( const DG1Function< real_t >&          src,
               const P2PlusBubbleFunction< real_t >& dst,
               size_t                                level,
               DoFType                               flag,
               UpdateType                            updateType ) const override
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "DG1ToP2PlusBubbleOperator::apply() does not support 3D meshes, yet." );
      }

      assembleAndOrApply2D< false >( src, dst, level, flag, nullptr, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const DG1Function< idx_t >&                 src,
                  const P2PlusBubbleFunction< idx_t >&        dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "DG1ToP2PlusBubbleOperator::toMatrix() does not support 3D meshes, yet." );
      }

      assembleAndOrApply2D< true >( src, dst, level, flag, mat, Replace );
   }

 private:
   template < bool assembleMatrix, typename value_t >
   inline void assembleAndOrApply2D( const DG1Function< value_t >&               src,
                                     const P2PlusBubbleFunction< value_t >&      dst,
                                     size_t                                      level,
                                     DoFType                                     flag,
                                     const std::shared_ptr< SparseMatrixProxy >& mat,
                                     UpdateType                                  updateType = Replace ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells() == false );

      if constexpr ( assembleMatrix )
      {
         WALBERLA_ASSERT( mat != nullptr );
      }
      else
      {
         WALBERLA_ASSERT( mat == nullptr );
      }

      DGBasisLinearLagrange_Example dstBasis;

      using indexing::Index;
      using volumedofspace::indexing::ElementNeighborInfo;

      const auto storage = this->getStorage();

      int dim = 2;

      // check, if we need to update stuff for the DG1 source function
      if ( storage->getAdditionalHaloDepth() > 0 )
      {
         src.communicate( level );
      }

      if ( updateType == Replace && mat == nullptr )
      {
         // We need to zero the destination array (including halos).
         // However, we must not zero out anything that is not flagged with the specified BCs.
         // Therefore we first zero out everything that is flagged, and then, later,
         // the halos of the highest dim primitives.
         dst.interpolate( real_c( 0 ), level, flag );
      }

      for ( auto& it : storage_->getFaces() )
      {
         Face&             face   = *it.second;
         const PrimitiveID faceID = it.first;

         const auto srcPolyDegree = 1;

         const auto numSrcDofs = 3;
         const auto numDstDofs = 7;

         const auto srcDofMemory = src.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );
         const auto srcMemLayout = src.getDGFunction()->volumeDoFFunction()->memoryLayout();

         PrimitiveDataID< FunctionMemory< value_t >, Face > dstVertexDoFID     = dst.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< value_t >, Face > dstEdgeDoFID       = dst.getEdgeDoFFunction().getFaceDataID();
         const auto                                         dstBubbleMemLayout = dst.getVolumeDoFFunction().memoryLayout();

         value_t* dstVertexData = face.getData( dstVertexDoFID )->getPointer( level );
         value_t* dstEdgeData   = face.getData( dstEdgeDoFID )->getPointer( level );
         value_t* dstBubbleData = dst.getVolumeDoFFunction().dofMemory( faceID, level );

         // zero out halos for matrix-free application
         if constexpr ( !assembleMatrix )
         {
            for ( const auto& idx : vertexdof::macroface::Iterator( level ) )
            {
               if ( vertexdof::macroface::isVertexOnBoundary( level, idx ) )
               {
                  auto arrayIdx           = vertexdof::macroface::index( level, idx.x(), idx.y() );
                  dstVertexData[arrayIdx] = real_c( 0 );
               }
            }

            for ( const auto& idx : edgedof::macroface::Iterator( level ) )
            {
               for ( const auto& orientation : edgedof::faceLocalEdgeDoFOrientations )
               {
                  if ( !edgedof::macroface::isInnerEdgeDoF( level, idx, orientation ) )
                  {
                     auto arrayIdx         = edgedof::macroface::index( level, idx.x(), idx.y(), orientation );
                     dstEdgeData[arrayIdx] = real_c( 0 );
                  }
               }
            }
         }

         // blue and gray faces
         for ( uint_t microFaceType = 0; microFaceType < 2; microFaceType++ )
         {
            auto faceType = facedof::allFaceTypes[microFaceType];
            auto itFace   = facedof::macroface::Iterator( level, faceType ).begin();

            while ( itFace != itFace.end() )
            {
               Index microFace = *itFace;
               itFace++;

               // TODO: all these coord computations can be executed _once_ and then the coordinates can be incremented by h
               // TODO: blending -> see issue #293

               // This object does the heavy lifting of computing all required coordinates and normals.
               ElementNeighborInfo neighborInfo;

               neighborInfo = ElementNeighborInfo( microFace, faceType, level, src.getBoundaryCondition(), faceID, storage_ );

               // Obtain micro-face's vertex coordinates
               const auto&              elementVertexCoords = neighborInfo.elementVertexCoords();
               std::array< Point3D, 3 > coords;
               coords[0] = elementVertexCoords[0];
               coords[1] = elementVertexCoords[1];
               coords[2] = elementVertexCoords[2];

               // Let the form calculate the local element matrix
               Matrixr< numDstDofs, numSrcDofs > localMat;
               form_->integrateAll( coords, localMat );

               // obtain indices into data-buffers of destination FE function
               std::array< uint_t, 3 > vertexDoFIndicesDst;
               vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndicesDst );

               std::array< uint_t, 3 > edgeDoFIndicesDst;
               edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndicesDst );

               uint_t bubbleDoFIndexDst =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 0, 1, level, dstBubbleMemLayout );

               // obtain indices into data-buffers of source FE function
               std::array< uint_t, 3 > dg1DoFIndicesSrc;
               dg1DoFIndicesSrc[0] =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 0, 3, level, srcMemLayout );
               dg1DoFIndicesSrc[1] =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 1, 3, level, srcMemLayout );
               dg1DoFIndicesSrc[2] =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 2, 3, level, srcMemLayout );

               // operator application
               if constexpr ( !assembleMatrix )
               {
                  Eigen::Matrix< real_t, numDstDofs, 1 > elVecDst;
                  Eigen::Matrix< real_t, numSrcDofs, 1 > elVecSrc;
                  elVecDst.setZero();
                  elVecSrc.setZero();

                  // fill local source vector
                  elVecSrc[0] = srcDofMemory[dg1DoFIndicesSrc[0]];
                  elVecSrc[1] = srcDofMemory[dg1DoFIndicesSrc[1]];
                  elVecSrc[2] = srcDofMemory[dg1DoFIndicesSrc[2]];

                  // matrix-vector multiplication
                  elVecDst = localMat * elVecSrc;

                  // distribute local destination vector
                  dstVertexData[vertexDoFIndicesDst[0]] += elVecDst[0];
                  dstVertexData[vertexDoFIndicesDst[1]] += elVecDst[1];
                  dstVertexData[vertexDoFIndicesDst[2]] += elVecDst[2];

                  dstEdgeData[edgeDoFIndicesDst[0]] += elVecDst[3];
                  dstEdgeData[edgeDoFIndicesDst[1]] += elVecDst[4];
                  dstEdgeData[edgeDoFIndicesDst[2]] += elVecDst[5];

                  dstBubbleData[bubbleDoFIndexDst] += elVecDst[6];
               }

               // matrix assembly
               else
               {
                  // vertex dofs
                  for ( int dstIdx = 0; dstIdx < 3; ++dstIdx )
                  {
                     const auto globalRowIdx = dstVertexData[vertexDoFIndicesDst[dstIdx]];

                     for ( int srcIdx = 0; srcIdx < 3; ++srcIdx )
                     {
                        const auto globalColIdx = srcDofMemory[dg1DoFIndicesSrc[srcIdx]];

                        mat->addValue( globalRowIdx, globalColIdx, localMat( dstIdx, srcIdx ) );
                     }
                  }

                  // edge dofs
                  for ( int dstIdx = 3; dstIdx < 6; ++dstIdx )
                  {
                     const auto globalRowIdx = dstEdgeData[edgeDoFIndicesDst[dstIdx - 3]];

                     for ( int srcIdx = 0; srcIdx < 3; ++srcIdx )
                     {
                        const auto globalColIdx = srcDofMemory[dg1DoFIndicesSrc[srcIdx]];

                        mat->addValue( globalRowIdx, globalColIdx, localMat( dstIdx, srcIdx ) );
                     }
                  }

                  // bubble dof
                  const auto globalRowIdx = dstBubbleData[bubbleDoFIndexDst];

                  for ( int srcIdx = 0; srcIdx < 3; ++srcIdx )
                  {
                     const auto globalColIdx = srcDofMemory[dg1DoFIndicesSrc[srcIdx]];

                     mat->addValue( globalRowIdx, globalColIdx, localMat( 6, srcIdx ) );
                  }
               }
            }
         }
      }

      if constexpr ( !assembleMatrix )
      {
         dst.getVertexDoFFunction().template communicateAdditively< Face, Edge >(
             level, DoFType::All ^ flag, *storage_, updateType == Replace );
         dst.getVertexDoFFunction().template communicateAdditively< Face, Vertex >(
             level, DoFType::All ^ flag, *storage_, updateType == Replace );
         dst.getEdgeDoFFunction().template communicateAdditively< Face, Edge >(
             level, DoFType::All ^ flag, *storage_, updateType == Replace );

         // as we never access the bubble dofs in the faces' halos we also do not need to communicate them here
      }
   }

   std::shared_ptr< Form > form_;
};

// DG1 to P2PlusBubble "gradient" (transpose divergence)
typedef DG1ToP2PlusBubbleOperator< forms::dg1_to_p2_plus_bubble_divt_0_affine_q3 > DG1ToP2PlusBubbleDivTxOperator;
typedef DG1ToP2PlusBubbleOperator< forms::dg1_to_p2_plus_bubble_divt_1_affine_q3 > DG1ToP2PlusBubbleDivTyOperator;
typedef DG1ToP2PlusBubbleOperator< forms::dg1_to_p2_plus_bubble_divt_2_affine_q3 > DG1ToP2PlusBubbleDivTzOperator;

} // namespace hyteg
