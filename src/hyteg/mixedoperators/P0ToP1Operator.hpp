/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/forms/form_hyteg_generated/p0_to_p1/p0_to_p1_divt_affine_q0.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

using namespace dg;
using facedof::FaceType;
using indexing::Index;
using volumedofspace::indexing::VolumeDoFMemoryLayout;
using walberla::int_c;
using walberla::real_t;

template < typename Form >
class P0ToP1Operator : public Operator< P0Function< real_t >, P1Function< real_t > >
{
 public:
   P0ToP1Operator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator< P0Function< real_t >, P1Function< real_t > >( storage, minLevel, maxLevel )
   , form_( std::make_shared< Form >() )
   {}

   void apply( const P0Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      assembleAndOrApply( src, dst, level, flag, nullptr, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P0Function< idx_t >&                  src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      assembleAndOrApply( src, dst, level, flag, mat, Replace );
   }

 private:
   /// \brief This is similar to the implementation in the dg::DGOperator class.
   template < typename VType >
   inline void assembleAndOrApply( const P0Function< VType >&                  src,
                                   const P1Function< VType >&                  dst,
                                   size_t                                      level,
                                   DoFType                                     flag,
                                   const std::shared_ptr< SparseMatrixProxy >& mat,
                                   UpdateType                                  updateType = Replace ) const
   {
      using indexing::Index;

      WALBERLA_CHECK( updateType == Replace );

      auto srcDGF = src.getDGFunction();

      if ( this->getStorage()->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P0 to P1 apply not implemented in 3D." );
      }
      else
      {
         const int dim = 2;

         for ( const auto& faceIt : this->getStorage()->getFaces() )
         {
            const auto  faceId = faceIt.first;
            const auto& face   = *faceIt.second;

            const auto srcPolyDegree = 0;
            const auto dstPolyDegree = 1;

            const auto numSrcDofs = 1;
            const auto numDstDofs = 3;

            const auto srcDofMemory = srcDGF->volumeDoFFunction()->dofMemory( faceId, level );
            auto       dstDofMemory = face.getData( dst.getFaceDataID() )->getPointer( level );

            const auto srcMemLayout = srcDGF->volumeDoFFunction()->memoryLayout();

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
               {
                  // TODO: all these coord computations can be executed _once_ and then the coordinates can be incremented by h
                  // TODO: blending

                  // This object does the heavy lifting of computing all required coordinates and normals.
                  volumedofspace::indexing::ElementNeighborInfo neighborInfo(
                      elementIdx, faceType, level, src.getBoundaryCondition(), faceId, this->getStorage() );

                  // We only write to the DoFs in the current volume, let's prepare a temporary vector for that.
                  Eigen::Matrix< real_t, Eigen::Dynamic, 1 > dstDofs;
                  dstDofs.resize( numDstDofs, Eigen::NoChange_t::NoChange );
                  dstDofs.setZero();

                  /////////////////////////
                  // Volume contribution //
                  /////////////////////////

                  std::array< Point3D, 3 > coords;
                  for ( uint_t i = 0; i < neighborInfo.elementVertexCoords().size(); i++ )
                  {
                     coords[i][0] = neighborInfo.elementVertexCoords().at( i )( 0 );
                     coords[i][1] = neighborInfo.elementVertexCoords().at( i )( 1 );
                     coords[i][2] = neighborInfo.elementVertexCoords().at( i )( 2 );
                  }

                  Matrixr< 3, 1 > elMat;
                  form_->integrateAll( coords, elMat );

                  // P0 has only one DoF
                  const auto srcDoF = srcDofMemory[volumedofspace::indexing::index(
                      elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, srcMemLayout )];

                  // Micro-vertex indices
                  auto vertexDoFIndices = facedof::macroface::getMicroVerticesFromMicroFace( elementIdx, faceType );

                  if ( mat == nullptr )
                  {
                     // Matrix-vector multiplication.
                     // dstDofs += localMat * srcDofs;
                     WALBERLA_ABORT( "Not implemented." );
                  }
                  else
                  {
                     // Sparse assembly.
                     for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
                     {
                        const auto globalRowIdx = dstDofMemory[vertexdof::macroface::index(
                            level, vertexDoFIndices.at( dstDofIdx ).x(), vertexDoFIndices.at( dstDofIdx ).y() )];
                        const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                            elementIdx.x(), elementIdx.y(), faceType, 0, 1, level, srcMemLayout )];
                        mat->addValue( globalRowIdx, globalColIdx, elMat( dstDofIdx, 0 ) );
                     }
                  }

                  if ( mat == nullptr )
                  {
                     // Write DoFs.
                     for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
                     {
                        WALBERLA_ABORT( "Not implemented." );
                     }
                  }
               }
            }
         }
      }

      WALBERLA_UNUSED( flag );
   }

   std::shared_ptr< Form > form_;
};

typedef P0ToP1Operator< forms::p0_to_p1_divt_0_affine_q0 > P0ToP1ConstantDivTxOperator;
typedef P0ToP1Operator< forms::p0_to_p1_divt_1_affine_q0 > P0ToP1ConstantDivTyOperator;
typedef P0ToP1Operator< forms::p0_to_p1_divt_2_affine_q0 > P0ToP1ConstantDivTzOperator;

} // namespace hyteg