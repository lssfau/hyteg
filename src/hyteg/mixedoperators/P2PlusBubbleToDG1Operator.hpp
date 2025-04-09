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
#include "hyteg/forms/form_hyteg_generated/p2_plus_bubble_to_dg1/p2_plus_bubble_to_dg1_div_affine_q3.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
// #include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2PlusBubbleFunction.hpp"

namespace hyteg {

using namespace dg;
using facedof::FaceType;
using indexing::Index;
using volumedofspace::indexing::VolumeDoFMemoryLayout;
using walberla::int_c;
using walberla::real_t;

template < typename Form >
class P2PlusBubbleToDG1Operator : public Operator< P2PlusBubbleFunction< real_t >, DG1Function< real_t > >
{
 public:
   P2PlusBubbleToDG1Operator( const std::shared_ptr< PrimitiveStorage >& storage,
                              uint_t                                     minLevel,
                              uint_t                                     maxLevel,
                              std::shared_ptr< Form >                    form = std::make_shared< Form >() )
   : Operator< P2PlusBubbleFunction< real_t >, DG1Function< real_t > >( storage, minLevel, maxLevel )
   , form_( form )
   {}

   typedef Form FormType;

   void apply( const P2PlusBubbleFunction< real_t >& src,
               const DG1Function< real_t >&          dst,
               size_t                                level,
               DoFType                               flag,
               UpdateType                            updateType ) const override
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2PlusBubbleToDG1Operator::apply() does not support 3D meshes, yet." );
      }

      assembleAndOrApply2D< false >( src, dst, level, flag, nullptr, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2PlusBubbleFunction< idx_t >&        src,
                  const DG1Function< idx_t >&                 dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2PlusBubbleToDG1Operator::toMatrix() does not support 3D meshes, yet." );
      }

      assembleAndOrApply2D< true >( src, dst, level, flag, mat, Replace );
   }

 private:
   template < bool assembleMatrix, typename value_t >
   inline void assembleAndOrApply2D( const P2PlusBubbleFunction< value_t >&      src,
                                     const DG1Function< value_t >&               dst,
                                     size_t                                      level,
                                     DoFType                                     flag,
                                     const std::shared_ptr< SparseMatrixProxy >& mat,
                                     UpdateType                                  updateType = Replace ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells() == false );
      WALBERLA_ASSERT( updateType == Add || updateType == Replace );

      if constexpr ( assembleMatrix )
      {
         WALBERLA_ASSERT( mat != nullptr );
      }
      else
      {
         WALBERLA_ASSERT( mat == nullptr );
      }

      using indexing::Index;
      using volumedofspace::indexing::ElementNeighborInfo;

      const auto storage = this->getStorage();

      int dim = 2;

      std::vector< PrimitiveID > pIDs = storage->getFaceIDs();

      communication::syncFunctionBetweenPrimitives( src, level );

      for ( auto& it : storage_->getFaces() )
      {
         Face&             face   = *it.second;
         const PrimitiveID faceID = it.first;

         const auto dstPolyDegree = 1;

         const auto numSrcDofs = 7;
         const auto numDstDofs = 3;

         const auto dstDofMemory = dst.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );
         const auto dstMemLayout = dst.getDGFunction()->volumeDoFFunction()->memoryLayout();

         PrimitiveDataID< FunctionMemory< value_t >, Face > srcVertexDoFID     = src.getVertexDoFFunction().getFaceDataID();
         PrimitiveDataID< FunctionMemory< value_t >, Face > srcEdgeDoFID       = src.getEdgeDoFFunction().getFaceDataID();
         const auto                                         srcBubbleMemLayout = src.getVolumeDoFFunction().memoryLayout();

         value_t* srcVertexData = face.getData( srcVertexDoFID )->getPointer( level );
         value_t* srcEdgeData   = face.getData( srcEdgeDoFID )->getPointer( level );
         value_t* srcBubbleData = src.getVolumeDoFFunction().dofMemory( faceID, level );

         // blue and gray faces
         for ( uint_t microFaceType = 0; microFaceType < 2; microFaceType++ )
         {
            auto faceType = facedof::allFaceTypes[microFaceType];

            for ( auto itFace = facedof::macroface::Iterator( level, faceType ).begin(); itFace != itFace.end(); ++itFace )
            {
               Index microFace = *itFace;

               // TODO: all these coord computations can be executed _once_ and then the coordinates can be incremented by h
               // TODO: blending

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
               std::array< uint_t, 3 > dg1DoFIndicesDst;
               dg1DoFIndicesDst[0] =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 0, 3, level, dstMemLayout );
               dg1DoFIndicesDst[1] =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 1, 3, level, dstMemLayout );
               dg1DoFIndicesDst[2] =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 2, 3, level, dstMemLayout );

               // obtain indices into data-buffers of source FE function
               std::array< uint_t, 3 > vertexDoFIndicesSrc;
               vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndicesSrc );

               std::array< uint_t, 3 > edgeDoFIndicesSrc;
               edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndicesSrc );

               uint_t bubbleDoFIndexSrc =
                   volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 0, 1, level, srcBubbleMemLayout );

               // operator application
               if constexpr ( !assembleMatrix )
               {
                  Eigen::Matrix< real_t, numDstDofs, 1 > elVecDst;
                  Eigen::Matrix< real_t, numSrcDofs, 1 > elVecSrc;
                  elVecDst.setZero();
                  elVecSrc.setZero();

                  // fill local source vector
                  elVecSrc[0] = srcVertexData[vertexDoFIndicesSrc[0]];
                  elVecSrc[1] = srcVertexData[vertexDoFIndicesSrc[1]];
                  elVecSrc[2] = srcVertexData[vertexDoFIndicesSrc[2]];

                  elVecSrc[3] = srcEdgeData[edgeDoFIndicesSrc[0]];
                  elVecSrc[4] = srcEdgeData[edgeDoFIndicesSrc[1]];
                  elVecSrc[5] = srcEdgeData[edgeDoFIndicesSrc[2]];

                  elVecSrc[6] = srcBubbleData[bubbleDoFIndexSrc];

                  // matrix-vector multiplication
                  elVecDst = localMat * elVecSrc;

                  // distribute local destination vector
                  if ( updateType == Replace )
                  {
                     dstDofMemory[dg1DoFIndicesDst[0]] = elVecDst[0];
                     dstDofMemory[dg1DoFIndicesDst[1]] = elVecDst[1];
                     dstDofMemory[dg1DoFIndicesDst[2]] = elVecDst[2];
                  }
                  else
                  {
                     dstDofMemory[dg1DoFIndicesDst[0]] += elVecDst[0];
                     dstDofMemory[dg1DoFIndicesDst[1]] += elVecDst[1];
                     dstDofMemory[dg1DoFIndicesDst[2]] += elVecDst[2];
                  }
               }

               // matrix assembly
               else
               {
                  // vertex dofs
                  for ( int srcIdx = 0; srcIdx < 3; ++srcIdx )
                  {
                     const auto globalColIdx = srcVertexData[vertexDoFIndicesSrc[srcIdx]];

                     for ( int dstIdx = 0; dstIdx < 3; ++dstIdx )
                     {
                        const auto globalRowIdx = dstDofMemory[dg1DoFIndicesDst[dstIdx]];

                        mat->addValue( globalRowIdx, globalColIdx, localMat( dstIdx, srcIdx ) );
                     }
                  }

                  // edge dofs
                  for ( int srcIdx = 3; srcIdx < 6; ++srcIdx )
                  {
                     const auto globalColIdx = srcEdgeData[edgeDoFIndicesSrc[srcIdx - 3]];

                     for ( int dstIdx = 0; dstIdx < 3; ++dstIdx )
                     {
                        const auto globalRowIdx = dstDofMemory[dg1DoFIndicesDst[dstIdx]];

                        mat->addValue( globalRowIdx, globalColIdx, localMat( dstIdx, srcIdx ) );
                     }
                  }

                  // bubble dof
                  const auto globalColIdx = srcBubbleData[bubbleDoFIndexSrc];

                  for ( int dstIdx = 0; dstIdx < 3; ++dstIdx )
                  {
                     const auto globalRowIdx = dstDofMemory[dg1DoFIndicesDst[dstIdx]];

                     mat->addValue( globalRowIdx, globalColIdx, localMat( dstIdx, 6 ) );
                  }
               }
            }
         }
      }
   }

   std::shared_ptr< Form > form_;
};

// P2PlusBubble to DG1 Stokes divergence
typedef P2PlusBubbleToDG1Operator< forms::p2_plus_bubble_to_dg1_div_0_affine_q3 > P2PlusBubbleToDG1DivxOperator;
typedef P2PlusBubbleToDG1Operator< forms::p2_plus_bubble_to_dg1_div_1_affine_q3 > P2PlusBubbleToDG1DivyOperator;
typedef P2PlusBubbleToDG1Operator< forms::p2_plus_bubble_to_dg1_div_2_affine_q3 > P2PlusBubbleToDG1DivzOperator;

} // namespace hyteg
