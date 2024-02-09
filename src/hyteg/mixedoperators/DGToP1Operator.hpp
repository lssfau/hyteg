/*
* Copyright (c) 2022 Andreas Wagner.
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
#include "hyteg/dgfunctionspace/DGFormAbort.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/P1_to_P0_div_form.hpp"
#include "hyteg/egfunctionspace/EGConstEpsilonForm.hpp"
#include "hyteg/egfunctionspace/EGDivForm.hpp"
#include "hyteg/egfunctionspace/EGDivtForm.hpp"
#include "hyteg/egfunctionspace/EGMassForm.hpp"
#include "hyteg/egfunctionspace/EGVectorLaplaceForm.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/mixedoperators/P1ToDGOperator.hpp"
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

template < typename Form, typename ValueType = real_t >
class DGToP1Operator : public Operator< DGFunction< ValueType >, P1Function< ValueType > >
{
 public:
   DGToP1Operator( const std::shared_ptr< PrimitiveStorage >& storage,
                   uint_t                                     minLevel,
                   uint_t                                     maxLevel,
                   std::shared_ptr< Form >                    form = std::make_shared< Form >() )
   : Operator< DGFunction< ValueType >, P1Function< ValueType > >( storage, minLevel, maxLevel )
   , form_( form )
   {}

   typedef Form FormType;

   void apply( const DGFunction< ValueType >& src,
               const P1Function< ValueType >& dst,
               size_t                         level,
               DoFType                        flag,
               UpdateType                     updateType ) const override
   {
      assembleAndOrApply( src, dst, level, flag, nullptr, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const DGFunction< idx_t >&                  src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      assembleAndOrApply( src, dst, level, flag, mat, Replace );
   }

 private:
   template < typename VType >
   VType* p1Data( const P1Function< VType >&                 function,
                  const std::shared_ptr< PrimitiveStorage >& storage,
                  const PrimitiveID&                         pid,
                  uint_t                                     level ) const
   {
      if ( storage->hasGlobalCells() )
      {
         WALBERLA_ASSERT( storage->cellExistsLocally( pid ) );
         auto cell = storage->getCell( pid );
         return cell->getData( function.getCellDataID() )->getPointer( level );
      }
      else
      {
         WALBERLA_ASSERT( storage->faceExistsLocally( pid ) );
         auto face = storage->getFace( pid );
         return face->getData( function.getFaceDataID() )->getPointer( level );
      }
   }

   /// \brief This is similar to the implementation in the dg::DGOperator class.
   template < typename VType >
   inline void assembleAndOrApply( const DGFunction< VType >&                  src,
                                   const P1Function< VType >&                  dst,
                                   size_t                                      level,
                                   DoFType                                     flag,
                                   const std::shared_ptr< SparseMatrixProxy >& mat,
                                   UpdateType                                  updateType ) const
   {
      // To avoid code duplication in this already long method, the implementation "fuses" the 2D and 3D implementation.
      // This more or less serves as a reference - for better performance the matrix-vector multiplication should be specialized.

      // zero out destination

      //VTKOutput vtk( "../../output", "DGToP1OpDebug", dst.getStorage() );
      //vtk.add( dst );
      //vtk.write( level, 0 );
      // vtk.write( level, 1 );

    //  WALBERLA_ASSERT(updateType == Replace);

      if (updateType == Replace && mat == nullptr)
         dst.interpolate(0, level, All);


      DGBasisLinearLagrange_Example dstBasis;

      using indexing::Index;
      using volumedofspace::indexing::ElementNeighborInfo;
      //communication::syncFunctionBetweenPrimitives( dst, level );
      src.communicate( level );

       const auto storage = this->getStorage();
       const int dim = storage->hasGlobalCells() ? 3 : 2;



      const std::vector< PrimitiveID > pids{ ( dim == 2 ) ? storage->getFaceIDs() : storage->getCellIDs() };

      for ( const auto& pid : pids )
      {
         const auto dstPolyDegree = 1;
         const auto srcPolyDegree = src.polynomialDegree( pid );

         const auto numDstDofs = dim + 1;
         const auto numSrcDofs = src.basis()->numDoFsPerElement( dim, dstPolyDegree );

         const auto dstDofMemory = p1Data< VType >( dst, storage, pid, level );
         auto       srcDofMemory = src.volumeDoFFunction()->dofMemory( pid, level );

         const auto srcMemLayout = src.volumeDoFFunction()->memoryLayout();

         // zero out halos for matrix-free application
         if ( mat == nullptr)
         {
            if ( dim == 2 )
            {
               for ( const auto& idx : vertexdof::macroface::Iterator( level ) )
               {
                  if ( vertexdof::macroface::isVertexOnBoundary( level, idx ) )
                  {
                     auto arrayIdx          = vertexdof::macroface::index( level, idx.x(), idx.y() );
                     dstDofMemory[arrayIdx] = real_c( 0 );
                  }
               }
            }
            else
            {
               for ( const auto& idx : vertexdof::macrocell::Iterator( level ) )
               {
                  if ( !vertexdof::macrocell::isOnCellFace( idx, level ).empty() )
                  {
                     auto arrayIdx          = vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() );
                     dstDofMemory[arrayIdx] = real_c( 0 );
                  }
               }
            }
         }

         const uint_t numMicroVolTypes = ( storage->hasGlobalCells() ? 6 : 2 );

         for ( uint_t microVolType = 0; microVolType < numMicroVolTypes; microVolType++ )
         {
            if ( dim == 2 && microVolType >= 2 )
            {
               break;
            }

            auto faceType = facedof::allFaceTypes[microVolType];
            auto cellType = celldof::allCellTypes[microVolType];

            auto itFace = facedof::macroface::Iterator( level, faceType ).begin();
            auto itCell = celldof::macrocell::Iterator( level, cellType ).begin();

            while ( ( dim == 2 && itFace != itFace.end() ) || ( dim == 3 && itCell != itCell.end() ) )
            {
               Index elementIdx;

               if ( dim == 2 )
               {
                  elementIdx = *itFace;
                  itFace++;
               }
               else
               {
                  elementIdx = *itCell;
                  itCell++;
               }

               // This object does the heavy lifting of computing all required coordinates and normals.
               ElementNeighborInfo neighborInfo;

               if ( dim == 2 )
               {
                  neighborInfo = ElementNeighborInfo( elementIdx, faceType, level, dst.getBoundaryCondition(), pid, storage );
               }
               else
               {
                  neighborInfo = ElementNeighborInfo( elementIdx, cellType, level, dst.getBoundaryCondition(), pid, storage );
               }

               // We only write to the DoFs in the current volume, let's prepare a temporary vector for that.
               PointXr dstDofs;
               dstDofs.resize( numDstDofs, Eigen::NoChange_t::NoChange );
               dstDofs.setZero();

               /////////////////////////
               // Volume contribution //
               /////////////////////////

               MatrixXr localMat;
               localMat.resize( numDstDofs, numSrcDofs );

               // Little difference here is that the source is now a CG P1 function.
               // So we need to obtain the DoFs a little differently and set the basis manually.

               form_->integrateVolume(
                   dim, neighborInfo.elementVertexCoords(), *src.basis(), dstBasis, srcPolyDegree, dstPolyDegree, localMat );

               PointXr srcDofs;
               srcDofs.resize( numSrcDofs, Eigen::NoChange_t::NoChange );

               // Getting the vertex DoF indices for the current micro volume.
               std::vector< Index > vertexDoFIndices;
               if ( dim == 2 )
               {
                  auto vertexDoFIndicesArray = facedof::macroface::getMicroVerticesFromMicroFace( elementIdx, faceType );
                  vertexDoFIndices.insert( vertexDoFIndices.begin(), vertexDoFIndicesArray.begin(), vertexDoFIndicesArray.end() );
               }
               else
               {
                  auto vertexDoFIndicesArray = celldof::macrocell::getMicroVerticesFromMicroCell( elementIdx, cellType );
                  vertexDoFIndices.insert( vertexDoFIndices.begin(), vertexDoFIndicesArray.begin(), vertexDoFIndicesArray.end() );
               }


               for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
               {
                  if ( dim == 2 )
                  {
                     srcDofs( srcDofIdx ) = srcDofMemory[volumedofspace::indexing::index(
                         elementIdx.x(), elementIdx.y(), faceType, srcDofIdx, numSrcDofs, level, srcMemLayout )];
                  }
                  else
                  {
                     srcDofs( srcDofIdx ) = srcDofMemory[volumedofspace::indexing::index(
                         elementIdx.x(), elementIdx.y(), elementIdx.z(), cellType, srcDofIdx, numSrcDofs, level, srcMemLayout )];
                  }
               }

               if ( mat == nullptr )
               {
                  // Matrix-vector multiplication.
                  dstDofs += localMat * srcDofs;
               }
               else
               {
                  // Sparse assembly.
                  for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
                  {
                     for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
                     {
                        if ( dim == 2 )
                        {
                           const auto globalRowIdx = dstDofMemory[vertexdof::macroface::index(
                               level, vertexDoFIndices[dstDofIdx].x(), vertexDoFIndices[dstDofIdx].y() )];
                           const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                               elementIdx.x(), elementIdx.y(), faceType, srcDofIdx, numSrcDofs, level, srcMemLayout )];
                           mat->addValue(
                               globalRowIdx, globalColIdx, localMat( Eigen::Index( dstDofIdx ), Eigen::Index( srcDofIdx ) ) );
                        }
                        else
                        {
                           const auto globalRowIdx = dstDofMemory[vertexdof::macrocell::index( level,
                                                                                               vertexDoFIndices[dstDofIdx].x(),
                                                                                               vertexDoFIndices[dstDofIdx].y(),
                                                                                               vertexDoFIndices[dstDofIdx].z() )];
                           const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index( elementIdx.x(),
                                                                                                   elementIdx.y(),
                                                                                                   elementIdx.z(),
                                                                                                   cellType,
                                                                                                   srcDofIdx,
                                                                                                   numSrcDofs,
                                                                                                   level,
                                                                                                   srcMemLayout )];
                           mat->addValue(
                               globalRowIdx, globalColIdx, localMat( Eigen::Index( dstDofIdx ), Eigen::Index( srcDofIdx ) ) );
                        }
                     }
                  }
               }

               if ( mat == nullptr )
               {
                  // Write DoFs.
                  for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
                  {
                     if ( dim == 2 )
                     {
                        const auto dstIdx = Eigen::Index( vertexdof::macroface::index(
                            level, vertexDoFIndices[dstDofIdx].x(), vertexDoFIndices[dstDofIdx].y() ) );

                        dstDofMemory[dstIdx] += dstDofs( dstDofIdx );
                     }
                     else
                     {
                        const auto dstIdx = Eigen::Index( vertexdof::macrocell::index( level,
                                                                                       vertexDoFIndices[dstDofIdx].x(),
                                                                                       vertexDoFIndices[dstDofIdx].y(),
                                                                                       vertexDoFIndices[dstDofIdx].z() ) );

                        dstDofMemory[dstIdx] += dstDofs( dstDofIdx );
                     }
                  }
               }
            }
         }


         if ( mat == nullptr )
         {
            if ( dim == 2 )
            {
               dst.template communicateAdditively< Face, Edge >( level, DoFType::All ^ flag, *storage, updateType == Replace );
               dst.template communicateAdditively< Face, Vertex >( level, DoFType::All ^ flag, *storage, updateType == Replace );
            }
            else
            {
               dst.template communicateAdditively< Cell, Face >( level, DoFType::All ^ flag, *storage, updateType == Replace );
               dst.template communicateAdditively< Cell, Edge >( level, DoFType::All ^ flag, *storage, updateType == Replace );
               dst.template communicateAdditively< Cell, Vertex >( level, DoFType::All ^ flag, *storage, updateType == Replace );
            }
         }
      }


   }

   std::shared_ptr< Form > form_;
};

} // namespace hyteg