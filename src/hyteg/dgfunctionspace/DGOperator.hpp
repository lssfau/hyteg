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
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {
namespace dg {

using celldof::CellType;
using facedof::FaceType;
using indexing::Index;
using volumedofspace::indexing::VolumeDoFMemoryLayout;
using walberla::int_c;
using walberla::real_t;

class DGOperator : public Operator< DGFunction< real_t >, DGFunction< real_t > >
{
 public:
   DGOperator( const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               const std::shared_ptr< DGForm >&           form );

   void apply( const DGFunction< real_t >& src,
               const DGFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const DGFunction< idx_t >&                  src,
                  const DGFunction< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override;

 private:
   /// Just a small helper method that writes the local matrix into the global sparse system.
   template < typename VType >
   void addLocalToGlobalMatrix( int                                                            dim,
                                int                                                            numSrcDofs,
                                int                                                            numDstDofs,
                                VType*                                                         srcDofMemory,
                                VType*                                                         dstDofMemory,
                                VolumeDoFMemoryLayout                                          srcMemLayout,
                                VolumeDoFMemoryLayout                                          dstMemLayout,
                                Index                                                          srcElementIdx,
                                Index                                                          dstElementIdx,
                                uint_t                                                         srcMicroVolType,
                                uint_t                                                         dstMicroVolType,
                                uint_t                                                         level,
                                std::shared_ptr< SparseMatrixProxy >                           mat,
                                const Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const
   {
      // Sparse assembly.
      if ( dim == 2 )
      {
         WALBERLA_ASSERT_LESS( srcMicroVolType, 2 );
         WALBERLA_ASSERT_LESS( dstMicroVolType, 2 );
         auto srcFaceType = facedof::allFaceTypes.at( srcMicroVolType );
         auto dstFaceType = facedof::allFaceTypes.at( dstMicroVolType );
         for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
         {
            for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
            {
               const auto globalRowIdx = dstDofMemory[volumedofspace::indexing::index(
                   dstElementIdx.x(), dstElementIdx.y(), dstFaceType, dstDofIdx, numDstDofs, level, dstMemLayout )];
               const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index(
                   srcElementIdx.x(), srcElementIdx.y(), srcFaceType, srcDofIdx, numSrcDofs, level, srcMemLayout )];
               mat->addValue( globalRowIdx, globalColIdx, localMat( dstDofIdx, srcDofIdx ) );
            }
         }
      }
      else
      {
         auto srcCellType = celldof::allCellTypes[srcMicroVolType];
         auto dstCellType = celldof::allCellTypes[dstMicroVolType];
         for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
         {
            for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
            {
               const auto globalRowIdx = dstDofMemory[volumedofspace::indexing::index( dstElementIdx.x(),
                                                                                       dstElementIdx.y(),
                                                                                       dstElementIdx.z(),
                                                                                       dstCellType,
                                                                                       dstDofIdx,
                                                                                       numDstDofs,
                                                                                       level,
                                                                                       dstMemLayout )];
               const auto globalColIdx = srcDofMemory[volumedofspace::indexing::index( srcElementIdx.x(),
                                                                                       srcElementIdx.y(),
                                                                                       srcElementIdx.z(),
                                                                                       srcCellType,
                                                                                       srcDofIdx,
                                                                                       numSrcDofs,
                                                                                       level,
                                                                                       srcMemLayout )];
               mat->addValue( globalRowIdx, globalColIdx, localMat( dstDofIdx, srcDofIdx ) );
            }
         }
      }
   }

   /// \brief Helper function that can be used to either apply the operator, or to assemble the sparse matrix.
   ///
   /// Since about 99% of the implementation is equal, it makes sense to fuse that here.
   ///
   /// Which operation is performed depends on the pointer to the sparse matrix proxy and the value type.
   /// If the pointer is a nullptr, apply is executed.
   ///
   /// Some notes on the implementation:
   ///
   /// For DG implementations there are two main possibilities to evaluate the interface integrals.
   ///
   /// The "naive" (not necessarily worse) approach is to loop over all interfaces in a dedicated loop e.g. after evaluating all
   /// volume integrals. Alternatively, the interface integrals are evaluated during the loop over the volumes.
   ///
   /// This function implements the latter, with the advantage that each DoF is only written to exactly once, also there is only
   /// a single loop over the macro-volume. On the downside, each interface integral has to be evaluated twice.
   ///
   /// A nice description is found in
   ///
   /// Kronbichler, M., & Kormann, K. (2019). Fast matrix-free evaluation of discontinuous Galerkin finite element operators.
   /// ACM Transactions on Mathematical Software (TOMS), 45(3), 1-40.
   ///
   template < typename VType >
   inline void assembleAndOrApply( const DGFunction< VType >&                  src,
                                   const DGFunction< VType >&                  dst,
                                   size_t                                      level,
                                   DoFType                                     flag,
                                   const std::shared_ptr< SparseMatrixProxy >& mat,
                                   UpdateType                                  updateType = Replace ) const
   {
      // To avoid code duplication in this already long method, the implementation "fuses" the 2D and 3D implementation.
      // This more or less serves as a reference - for better performance the matrix-vector multiplication should be specialized.

      using indexing::Index;
      using volumedofspace::indexing::ElementNeighborInfo;

      WALBERLA_CHECK( updateType == Replace );

      if ( !form_->onlyVolumeIntegrals() )
      {
         src.communicate( level );
      }

      const auto storage = this->getStorage();

      int dim = 2;
      if ( storage->hasGlobalCells() )
      {
         dim = 3;
      }

      std::vector< PrimitiveID > pids;
      if ( dim == 2 )
      {
         pids = storage->getFaceIDs();
      }
      else
      {
         pids = storage->getCellIDs();
      }

      for ( const auto& pid : pids )
      {
         const auto srcPolyDegree = src.polynomialDegree( pid );
         const auto dstPolyDegree = dst.polynomialDegree( pid );

         const auto numSrcDofs = src.basis()->numDoFsPerElement( dim, srcPolyDegree );
         const auto numDstDofs = dst.basis()->numDoFsPerElement( dim, dstPolyDegree );

         const auto srcDofMemory = src.volumeDoFFunction()->dofMemory( pid, level );
         auto       dstDofMemory = dst.volumeDoFFunction()->dofMemory( pid, level );

         const auto srcMemLayout = src.volumeDoFFunction()->memoryLayout();
         const auto dstMemLayout = dst.volumeDoFFunction()->memoryLayout();

         std::map< uint_t, VType* > glMemory;

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

               // TODO: all these coord computations can be executed _once_ and then the coordinates can be incremented by h
               // TODO: blending

               // This object does the heavy lifting of computing all required coordinates and normals.
               ElementNeighborInfo neighborInfo;

               if ( dim == 2 )
               {
                  neighborInfo = ElementNeighborInfo( elementIdx, faceType, level, src.getBoundaryCondition(), pid, storage_ );
               }
               else
               {
                  neighborInfo = ElementNeighborInfo( elementIdx, cellType, level, src.getBoundaryCondition(), pid, storage_ );
               }

               // We only write to the DoFs in the current volume, let's prepare a temporary vector for that.
               Eigen::Matrix< real_t, Eigen::Dynamic, 1 > dstDofs;
               dstDofs.resize( numDstDofs, Eigen::NoChange_t::NoChange );
               dstDofs.setZero();

               /////////////////////////
               // Volume contribution //
               /////////////////////////

               Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > localMat;
               localMat.resize( numDstDofs, numSrcDofs );

               form_->integrateVolume(
                   dim, neighborInfo.elementVertexCoords(), *src.basis(), *dst.basis(), srcPolyDegree, dstPolyDegree, localMat );

               // Volume DoFs are source.
               Eigen::Matrix< real_t, Eigen::Dynamic, 1 > srcDofs;
               srcDofs.resize( numSrcDofs, Eigen::NoChange_t::NoChange );

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
                  addLocalToGlobalMatrix( dim,
                                          numSrcDofs,
                                          numDstDofs,
                                          srcDofMemory,
                                          dstDofMemory,
                                          srcMemLayout,
                                          dstMemLayout,
                                          elementIdx,
                                          elementIdx,
                                          microVolType,
                                          microVolType,
                                          level,
                                          mat,
                                          localMat );
               }

               if ( !form_->onlyVolumeIntegrals() )
               {
                  /////////////////////////////
                  // Interface contributions //
                  /////////////////////////////

                  // Loop over neighboring volumes.
                  for ( uint_t n = 0; n < uint_c( dim + 1 ); n++ )
                  {
                     /////////////////////
                     // Domain boundary //
                     /////////////////////

                     if ( neighborInfo.atMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == DirichletBoundary )
                     {
                        ////////////////////////
                        // Dirichlet boundary //
                        ////////////////////////

                        localMat.setZero();
                        form_->integrateFacetDirichletBoundary( dim,
                                                                neighborInfo.elementVertexCoords(),
                                                                neighborInfo.interfaceVertexCoords( n ),
                                                                neighborInfo.oppositeVertexCoords( n ),
                                                                neighborInfo.outwardNormal( n ),
                                                                *src.basis(),
                                                                *dst.basis(),
                                                                srcPolyDegree,
                                                                dstPolyDegree,
                                                                localMat );

                        if ( mat == nullptr )
                        {
                           // Matrix-vector multiplication.
                           dstDofs += localMat * srcDofs;
                        }
                        else
                        {
                           // Sparse assembly.
                           addLocalToGlobalMatrix( dim,
                                                   numSrcDofs,
                                                   numDstDofs,
                                                   srcDofMemory,
                                                   dstDofMemory,
                                                   srcMemLayout,
                                                   dstMemLayout,
                                                   elementIdx,
                                                   elementIdx,
                                                   microVolType,
                                                   microVolType,
                                                   level,
                                                   mat,
                                                   localMat );
                        }
                     }
                     else if ( neighborInfo.atMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == NeumannBoundary )
                     {
                        WALBERLA_ABORT( "Neumann boundary handling not implemented." );
                     }
                     else if ( neighborInfo.atMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == FreeslipBoundary )
                     {
                        WALBERLA_ABORT( "Free-slip boundary handling not implemented." );
                     }

                     //////////////////
                     // Inner domain //
                     //////////////////

                     // This is the place where we need to handle hanging nodes. Especially during the outer element coupling.
                     //
                     // The main task is to get the micro-element vertices in the same order as they are local to the neighbor
                     // macro.
                     // Also, we need to get the actual coordinates of those vertices. And some more info like the type of
                     // micro-volume (orientation of the face or cell) and the ghost-layer micro-element indices.
                     //
                     // There are three cases to be considered:
                     //
                     //   a) equal level macro-macro interface
                     //   b) the current macro is finer (ghost-layer is coarser)
                     //   c) the current macro is coarser (ghost-layer is finer)
                     //
                     // Note that in case (c) there are multiple neighbor elements that we need to take into account
                     // (2 in 2D, 4 in 3D). To resolve the indices we chain some index transformations to get to the actual
                     // indices from the perspective of the neighboring macro.
                     //
                     // In case b) we first translate the micro-index to that index in the parent macro with a higher
                     // micro-refinement level. From that point we can find the neighboring index for that micro-refinement level
                     // on the neighboring macro via case a). Then we find the parent micro-index of that.
                     //
                     // In case c) we do something similar but the other way around: first we pretend to be on a finer
                     // micro-refinement level (there are now multiple possible micro-elements). We loop over those micros that
                     // are at the boundary and get their neighbors as in case a). Then we find the corresponding micro-indices
                     // local to the neighbor macro on a finer macro-refinement level.
                     //
                     // Detailed comments in the implementation below.

                     else
                     {
                        ///////////////////////////////////
                        // a) inner element contribution //
                        ///////////////////////////////////

                        // Although not obvious at this point we have to cover the AMR case also for the inner integrals.
                        //
                        // This applies only if the current macro is coarser than the neighboring macro(s).
                        // The reason is that the integral over the interface is split into 2 (2D) or 4 (3D) integrals -
                        // one for each neighbor element. Since we split the outer coupling we need to split the inner coupling
                        // as well.

                        // Let's first check the refinement level of the neighboring macro, if there is any.

                        bool hasNeighbor = neighborInfo.atMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == Inner;

                        uint_t macroLevel, neighborMacroLevel;
                        if ( hasNeighbor )
                        {
                           if ( dim == 2 )
                           {
                              const auto nMacroFaceIDs =
                                  storage->getFace( pid )->getIndirectTopLevelNeighborFaceIDsOverEdges().at(
                                      neighborInfo.macroBoundaryID( n ) );
                              macroLevel         = storage->getRefinementLevel( pid );
                              neighborMacroLevel = storage->getRefinementLevel( nMacroFaceIDs[0] );
                           }
                           else
                           {
                              // TODO: This must be updated for 3D.
                              const auto nMacroCellID = storage->getCell( pid )->getIndirectNeighborCellIDsOverFaces().at(
                                  neighborInfo.macroBoundaryID( n ) );
                              macroLevel         = storage->getRefinementLevel( pid );
                              neighborMacroLevel = storage->getRefinementLevel( nMacroCellID );
                              WALBERLA_CHECK_EQUAL( macroLevel, neighborMacroLevel );
                           }
                        }

                        // We collect the interface(s) to integrate in this list.
                        std::vector< std::vector< Eigen::Matrix< real_t, 3, 1 > > > interfaceVertexCoordsList;

                        if ( !hasNeighbor || macroLevel == neighborMacroLevel || macroLevel == neighborMacroLevel + 1 )
                        {
                           // a) Equal level or b) neighbor is coarser.
                           // Current macro is on the right.
                           /*
                                   *           *
                                  /|\         /|\
                                 / | \       / | \
                                /  |  \     /  |  \
                               *   *   *   *   *---*
                                \  |  /     \  |
                                 \ | /       \ |
                                  \|/         \|
                                   *           *
                           */

                           // Nothing special has to be done. We simply append the single interface we integrate to the list.

                           interfaceVertexCoordsList.push_back( neighborInfo.interfaceVertexCoords( n ) );
                        }
                        else if ( macroLevel + 1 == neighborMacroLevel )
                        {
                           // c) Neighbor is finer.
                           //    Current macro is on the right.
                           //
                           /*
                                   *
                                  /|\
                                 / | \
                                /  |  \
                               *---*   *
                                \  |  /
                                 \ | /
                                  \|/
                                   *
                           */

                           // We need to split the interface accordingly and perform 2 (2D) or 4 (3D) integrations.
                           // The corresponding volume is the same in all cases. Only the interface changes.

                           if ( dim == 2 )
                           {
                              std::vector< Index >    fineElementIndices;
                              std::vector< FaceType > fineFaceTypes;

                              volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                                  elementIdx, faceType, fineElementIndices, fineFaceTypes );

                              WALBERLA_CHECK_EQUAL( fineElementIndices.size(), 4 );
                              WALBERLA_CHECK_EQUAL( fineFaceTypes.size(), 4 );

                              for ( uint_t i = 0; i < 4; i++ )
                              {
                                 // This guy is the neighborhood info for refinement level = level + 1 on the current macro.
                                 ElementNeighborInfo fakeNeighborInfo( fineElementIndices[i],
                                                                       fineFaceTypes[i],
                                                                       level + 1,
                                                                       src.getBoundaryCondition(),
                                                                       pid,
                                                                       storage );

                                 if ( fakeNeighborInfo.atMacroBoundary( n ) )
                                 {
                                    WALBERLA_CHECK_EQUAL( fakeNeighborInfo.faceType(), facedof::FaceType::GRAY );

                                    interfaceVertexCoordsList.push_back( fakeNeighborInfo.interfaceVertexCoords( n ) );
                                 }
                              }

                              WALBERLA_CHECK_EQUAL( interfaceVertexCoordsList.size(), 2 );
                           }
                           else
                           {
                              WALBERLA_ABORT( "Not implemented for 3D." );
                           }
                        }
                        else
                        {
                           WALBERLA_ABORT( "Invalid mesh refinement - 2:1 balance probably not satisfied." );
                        }

                        for ( uint_t i = 0; i < interfaceVertexCoordsList.size(); i++ )
                        {
                           localMat.setZero();
                           form_->integrateFacetInner( dim,
                                                       neighborInfo.elementVertexCoords(),
                                                       interfaceVertexCoordsList[i],
                                                       neighborInfo.oppositeVertexCoords( n ),
                                                       neighborInfo.outwardNormal( n ),
                                                       *src.basis(),
                                                       *dst.basis(),
                                                       srcPolyDegree,
                                                       dstPolyDegree,
                                                       localMat );

                           if ( mat == nullptr )
                           {
                              // Matrix-vector multiplication.
                              dstDofs += localMat * srcDofs;
                           }
                           else
                           {
                              // Sparse assembly.
                              addLocalToGlobalMatrix( dim,
                                                      numSrcDofs,
                                                      numDstDofs,
                                                      srcDofMemory,
                                                      dstDofMemory,
                                                      srcMemLayout,
                                                      dstMemLayout,
                                                      elementIdx,
                                                      elementIdx,
                                                      microVolType,
                                                      microVolType,
                                                      level,
                                                      mat,
                                                      localMat );
                           }
                        }
                        ////////////////////////////////////////
                        // b) coupling to neighboring element //
                        ////////////////////////////////////////

                        if ( neighborInfo.atMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == Inner )
                        {
                           ////////////////////////////////////////////////
                           // i) micro-interface on macro-macro-boundary //
                           ////////////////////////////////////////////////

                           // Here comes the hard work to get AMR running.

                           // The following arrays should all have the same length - one element per interface integral.
                           // Eventually fo cases (a) and (b) the vectors should only have one element.
                           // Only for case (c) they are larger (since there are multiple element interface to integrate over).

                           // Vertex coords of the neighbor elements.
                           std::vector< std::vector< Eigen::Matrix< real_t, 3, 1 > > > neighborElementVertexCoordsList;

                           // Interface-opposite vertex of the neighboring element..
                           std::vector< Eigen::Matrix< real_t, 3, 1 > > neighborElementOppositeVertexCoordsList;

                           // We must choose the interface that is "smaller", i.e. the interface of the finer element.
                           // That is the interface we need to integrate over.
                           std::vector< std::vector< Eigen::Matrix< real_t, 3, 1 > > > interfaceVertexCoordsList;

                           // Finally we need the corresponding array indices in the ghost-layer.
                           // To get those, we store the corresponding neighboring inner micro-element indices,
                           // the refinement level and local element cell type so that we can plug the correct stuff into the
                           // indexing function later.
                           std::vector< Index >    glInnerMicroElementIdx;
                           std::vector< FaceType > glInnerMicroElementFaceType;
                           std::vector< CellType > glInnerMicroElementCellType;
                           std::vector< uint_t >   glInnerLevel;

                           if ( macroLevel == neighborMacroLevel )
                           {
                              // a) Equal level.
                              //
                              /*
                                   *
                                  /|\
                                 / | \
                                /  |  \
                               *   *   *
                                \  |  /
                                 \ | /
                                  \|/
                                   *
                               */
                              //
                              // All necessary information can be pulled straight out of the neighbor element info.

                              glMemory[n] = src.volumeDoFFunction()->glMemory( pid, level, neighborInfo.macroBoundaryID( n ) );

                              std::vector< Eigen::Matrix< real_t, 3, 1 > > neighborElementVertexCoords;
                              Eigen::Matrix< real_t, 3, 1 >                neighborElementOppositeVertexCoords;

                              neighborInfo = neighborInfo.updateForMacroBoundary( n );

                              neighborElementVertexCoordsList.push_back( neighborInfo.neighborElementVertexCoords( n ) );
                              neighborElementOppositeVertexCoordsList.push_back( neighborInfo.neighborOppositeVertexCoords( n ) );
                              interfaceVertexCoordsList.push_back( neighborInfo.interfaceVertexCoords( n ) );

                              glInnerMicroElementIdx.push_back( elementIdx );
                              glInnerMicroElementFaceType.push_back( faceType );
                              glInnerMicroElementCellType.push_back( cellType );
                              glInnerLevel.push_back( level );
                           }
                           else if ( macroLevel == neighborMacroLevel + 1 )
                           {
                              // b) Neighbor is coarser.
                              //    Current macro is on the right.
                              //
                              /*
                                   *
                                  /|\
                                 / | \
                                /  |  \
                               *   *---*
                                \  |
                                 \ |
                                  \|
                                   *
                               */
                              //
                              // We emulate the correct neighborhood information by pretending that the
                              // local macro-element is coarser than it actually is. We find the corresponding micro-element
                              // of the next coarser macro-element and stuff this data into the neighbor element info.
                              // We then obtain the neighbor micro-element info as if we were in case (a).
                              //
                              // The interface is taken from the local micro-element, though since it is the smaller one.
                              //
                              // Again: there is only _one_ neighbor micro-element we care about.

                              if ( dim == 2 )
                              {
                                 glMemory[n] =
                                     src.volumeDoFFunction()->glMemory( pid, level - 1, neighborInfo.macroBoundaryID( n ) );

                                 PrimitiveID coarsePID;
                                 Index       tmpCoarseElementIdx;
                                 FaceType    tmpCoarseFaceType;

                                 volumedofspace::indexing::getVolumeIdxOnCoarseMacro( *storage,
                                                                                      *storage->getFace( pid ),
                                                                                      level,
                                                                                      elementIdx,
                                                                                      faceType,
                                                                                      coarsePID,
                                                                                      tmpCoarseElementIdx,
                                                                                      tmpCoarseFaceType );

                                 // Now we have the micro-volume idx for the coarse macro on refinement level = level + 1.
                                 // But we need it on refinement level = level.

                                 Index    coarseElementIdx;
                                 FaceType coarseFaceType;

                                 volumedofspace::indexing::getCoarseMicroElementFromFineMicroElement(
                                     tmpCoarseElementIdx, tmpCoarseFaceType, coarseElementIdx, coarseFaceType );

                                 // The neighborhood info emulates now that we are on a coarser macro, but the same refinement
                                 // level. This is already what we need.

                                 ElementNeighborInfo fakeNeighborInfo(
                                     coarseElementIdx, coarseFaceType, level, src.getBoundaryCondition(), coarsePID, storage );

                                 // We can now retrieve all the infos as in case (a).

                                 // I think the following is only true for Bey's ordering in 2D, in 3D we need to find the correct
                                 // neighbor of the "fake" macro.
                                 WALBERLA_ASSERT( fakeNeighborInfo.atMacroBoundary( n ) );
                                 fakeNeighborInfo = fakeNeighborInfo.updateForMacroBoundary( n );

                                 neighborElementVertexCoordsList.push_back( fakeNeighborInfo.neighborElementVertexCoords( n ) );
                                 neighborElementOppositeVertexCoordsList.push_back(
                                     fakeNeighborInfo.neighborOppositeVertexCoords( n ) );

                                 // The current element is finer, so we need to integrate over its (smaller) interface.
                                 // We get that from the neighborInfo computed from the _actual_ micro-element.
                                 interfaceVertexCoordsList.push_back( neighborInfo.interfaceVertexCoords( n ) );

                                 Index    glElementIdx;
                                 FaceType glFaceType;
                                 CellType glCellType = CellType::WHITE_UP; // not relevant for 2D
                                 volumedofspace::indexing::getCoarseMicroElementFromFineMicroElement(
                                     elementIdx, faceType, glElementIdx, glFaceType );

                                 glInnerMicroElementIdx.push_back( glElementIdx );
                                 glInnerMicroElementFaceType.push_back( glFaceType );
                                 glInnerMicroElementCellType.push_back( glCellType );
                                 glInnerLevel.push_back( level - 1 );
                              }
                              else
                              {
                                 WALBERLA_ABORT( "Not implemented for 3D." );
                              }
                           }
                           else if ( macroLevel + 1 == neighborMacroLevel )
                           {
                              // c) Neighbor is finer.
                              //    Current macro is on the right.
                              //
                              /*
                                   *
                                  /|\
                                 / | \
                                /  |  \
                               *---*   *
                                \  |  /
                                 \ | /
                                  \|/
                                   *
                               */
                              //
                              // This is the most complicated case since
                              //  - there are multiple interfaces to integrate and
                              //  - we cannot emulate macro-refinement (applying the algorithm from case (b) doesn't work directly).
                              // Instead, we
                              //  1. pretend that we are on refinement level = level + 1
                              //  2. find all local micros that overlap with the current micro on refinement level = level + 1
                              //  3. construct the neighbor info for those micros
                              //  4. reject all micros that are not on the boundary
                              //  5. use the neighbor micros (on the neighbor macro) and find their equivalents on the finer
                              //     (neighbor) macros

                              if ( dim == 2 )
                              {
                                 glMemory[n] =
                                     src.volumeDoFFunction()->glMemory( pid, level + 1, neighborInfo.macroBoundaryID( n ) );

                                 std::vector< Index >    fineElementIndices;
                                 std::vector< FaceType > fineFaceTypes;

                                 volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                                     elementIdx, faceType, fineElementIndices, fineFaceTypes );

                                 WALBERLA_CHECK_EQUAL( fineElementIndices.size(), 4 );
                                 WALBERLA_CHECK_EQUAL( fineFaceTypes.size(), 4 );

                                 for ( uint_t i = 0; i < 4; i++ )
                                 {
                                    // This guy is the neighborhood info for refinement level = level + 1 on the current macro.
                                    ElementNeighborInfo fakeNeighborInfo( fineElementIndices[i],
                                                                          fineFaceTypes[i],
                                                                          level + 1,
                                                                          src.getBoundaryCondition(),
                                                                          pid,
                                                                          storage );

                                    if ( fakeNeighborInfo.atMacroBoundary( n ) )
                                    {
                                       WALBERLA_CHECK_EQUAL( fakeNeighborInfo.faceType(), facedof::FaceType::GRAY );

                                       // Getting neighborhood information for the equal-level neighbor macro on refinement
                                       // level = level + 1.
                                       fakeNeighborInfo = fakeNeighborInfo.updateForMacroBoundary( n );

                                       PrimitiveID fineNeighborPID;
                                       Index       fineNeighborElementIdx;
                                       FaceType    fineNeighborFaceType;

                                       auto localFace            = storage->getFace( pid );
                                       auto neighborCoarseFaceID = localFace->getIndirectNeighborFaceIDsOverEdges().at(
                                           fakeNeighborInfo.macroBoundaryID( n ) );
                                       auto neighborCoarseFace = storage->getFace( neighborCoarseFaceID );

                                       // Getting the _correct_ element on refinement level = level on the _correct_ finer
                                       // neighbor macro.

                                       volumedofspace::indexing::getVolumeIdxOnRefinedMacro(
                                           *neighborCoarseFace,
                                           level,
                                           fakeNeighborInfo.neighborElementIndices( n ),
                                           fakeNeighborInfo.neighborFaceType( n ),
                                           fineNeighborPID,
                                           fineNeighborElementIdx,
                                           fineNeighborFaceType );

                                       const auto& fineNeighborFace = *storage->getFace( fineNeighborPID );

                                       std::vector< Eigen::Matrix< real_t, 3, 1 > > neighborElementVertexCoords( 3 );
                                       Eigen::Matrix< real_t, 3, 1 >                neighborElementOppositeVertexCoords;

                                       auto microVertexIndices = facedof::macroface::getMicroVerticesFromMicroFace(
                                           fineNeighborElementIdx, fineNeighborFaceType );
                                       for ( uint_t ii = 0; ii < 3; ii++ )
                                       {
                                          auto tmp = vertexdof::macroface::coordinateFromIndex(
                                              level, fineNeighborFace, microVertexIndices[ii] );
                                          neighborElementVertexCoords[ii]( 0 ) = tmp[0];
                                          neighborElementVertexCoords[ii]( 1 ) = tmp[1];
                                          neighborElementVertexCoords[ii]( 2 ) = 0;
                                       }
                                       neighborElementVertexCoordsList.push_back( neighborElementVertexCoords );

                                       // Since the element is the same, we can use the information from the updated fake neighbor
                                       // info. Whew, otherwise we would need some more code...
                                       neighborElementOppositeVertexCoordsList.push_back(
                                           fakeNeighborInfo.oppositeVertexCoords( n ) );

                                       // The neighbor element is finer, so need to integrate over its interface.
                                       // The fake neighbor info is already refined, so we can use that interface.
                                       interfaceVertexCoordsList.push_back( fakeNeighborInfo.interfaceVertexCoords( n ) );

                                       glInnerMicroElementIdx.push_back( fineElementIndices[i] );
                                       glInnerMicroElementFaceType.push_back( fineFaceTypes[i] );
                                       glInnerMicroElementCellType.push_back( CellType::WHITE_UP ); // to be fixed for 3D
                                       glInnerLevel.push_back( level + 1 );
                                    }
                                 }

                                 WALBERLA_CHECK_EQUAL( interfaceVertexCoordsList.size(), 2 );
                              }
                              else
                              {
                                 WALBERLA_ABORT( "Not implemented for 3D." );
                              }
                           }
                           else
                           {
                              WALBERLA_ABORT( "Invalid mesh refinement - 2:1 balance probably not satisfied." );
                           }

                           WALBERLA_CHECK_EQUAL( neighborElementVertexCoordsList.size(),
                                                 neighborElementOppositeVertexCoordsList.size() );
                           WALBERLA_CHECK_EQUAL( neighborElementVertexCoordsList.size(), interfaceVertexCoordsList.size() );
                           WALBERLA_CHECK_EQUAL( neighborElementVertexCoordsList.size(), glInnerMicroElementIdx.size() );
                           WALBERLA_CHECK_EQUAL( neighborElementVertexCoordsList.size(), glInnerMicroElementFaceType.size() );
                           WALBERLA_CHECK_EQUAL( neighborElementVertexCoordsList.size(), glInnerMicroElementCellType.size() );
                           WALBERLA_CHECK_EQUAL( neighborElementVertexCoordsList.size(), glInnerLevel.size() );

                           // The info vectors only have one element in cases (a) and (b), but 2 or 4 in case (c).
                           for ( uint_t i = 0; i < neighborElementVertexCoordsList.size(); i++ )
                           {
                              const auto interfaceVertexCoords = interfaceVertexCoordsList[i];

                              localMat.setZero();
                              form_->integrateFacetCoupling( dim,
                                                             neighborInfo.elementVertexCoords(),
                                                             neighborElementVertexCoordsList[i],
                                                             interfaceVertexCoords,
                                                             neighborInfo.oppositeVertexCoords( n ),
                                                             neighborElementOppositeVertexCoordsList[i],
                                                             neighborInfo.outwardNormal( n ),
                                                             *src.basis(),
                                                             *dst.basis(),
                                                             srcPolyDegree,
                                                             dstPolyDegree,
                                                             localMat );

                              // Now we need the DoFs from the neighboring element.
                              Eigen::Matrix< real_t, Eigen::Dynamic, 1 > nSrcDofs;
                              nSrcDofs.resize( numSrcDofs, Eigen::NoChange_t::NoChange );
                              std::vector< uint_t > nSrcDoFArrIndices( numSrcDofs );

                              for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
                              {
                                 if ( dim == 2 )
                                 {
                                    nSrcDoFArrIndices[srcDofIdx] =
                                        volumedofspace::indexing::indexNeighborInGhostLayer( neighborInfo.macroBoundaryID( n ),
                                                                                             glInnerMicroElementIdx[i].x(),
                                                                                             glInnerMicroElementIdx[i].y(),
                                                                                             glInnerMicroElementFaceType[i],
                                                                                             srcDofIdx,
                                                                                             numSrcDofs,
                                                                                             glInnerLevel[i],
                                                                                             srcMemLayout );
                                 }
                                 else
                                 {
                                    nSrcDoFArrIndices[srcDofIdx] =
                                        volumedofspace::indexing::indexNeighborInGhostLayer( neighborInfo.macroBoundaryID( n ),
                                                                                             glInnerMicroElementIdx[i].x(),
                                                                                             glInnerMicroElementIdx[i].y(),
                                                                                             glInnerMicroElementIdx[i].z(),
                                                                                             glInnerMicroElementCellType[i],
                                                                                             srcDofIdx,
                                                                                             numSrcDofs,
                                                                                             glInnerLevel[i],
                                                                                             srcMemLayout );
                                 }

                                 nSrcDofs( srcDofIdx ) = glMemory[n][nSrcDoFArrIndices[srcDofIdx]];
                              }

                              if ( mat == nullptr )
                              {
                                 // Matrix-vector multiplication.
                                 dstDofs += localMat * nSrcDofs;
                              }
                              else
                              {
                                 // Sparse assembly.
                                 // TODO: maybe there is a nicer way to do the gl stuff ...
                                 for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
                                 {
                                    for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
                                    {
                                       uint_t globalRowIdx;
                                       if ( dim == 2 )
                                       {
                                          globalRowIdx = dstDofMemory[volumedofspace::indexing::index( elementIdx.x(),
                                                                                                       elementIdx.y(),
                                                                                                       faceType,
                                                                                                       dstDofIdx,
                                                                                                       numDstDofs,
                                                                                                       level,
                                                                                                       dstMemLayout )];
                                       }
                                       else
                                       {
                                          globalRowIdx = dstDofMemory[volumedofspace::indexing::index( elementIdx.x(),
                                                                                                       elementIdx.y(),
                                                                                                       elementIdx.z(),
                                                                                                       cellType,
                                                                                                       dstDofIdx,
                                                                                                       numDstDofs,
                                                                                                       level,
                                                                                                       dstMemLayout )];
                                       }
                                       const auto globalColIdx = glMemory[n][nSrcDoFArrIndices[srcDofIdx]];

                                       mat->addValue( globalRowIdx, globalColIdx, localMat( dstDofIdx, srcDofIdx ) );
                                    }
                                 }
                              }
                           }
                        }
                        else
                        {
                           /////////////////////////////////////////
                           // ii) micro-interface inside of macro //
                           /////////////////////////////////////////

                           localMat.setZero();
                           form_->integrateFacetCoupling( dim,
                                                          neighborInfo.elementVertexCoords(),
                                                          neighborInfo.neighborElementVertexCoords( n ),
                                                          neighborInfo.interfaceVertexCoords( n ),
                                                          neighborInfo.oppositeVertexCoords( n ),
                                                          neighborInfo.neighborOppositeVertexCoords( n ),
                                                          neighborInfo.outwardNormal( n ),
                                                          *src.basis(),
                                                          *dst.basis(),
                                                          srcPolyDegree,
                                                          dstPolyDegree,
                                                          localMat );

                           // Now we need the DoFs from the neighboring element.
                           Eigen::Matrix< real_t, Eigen::Dynamic, 1 > nSrcDofs;
                           nSrcDofs.resize( numSrcDofs, Eigen::NoChange_t::NoChange );

                           for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
                           {
                              if ( dim == 2 )
                              {
                                 nSrcDofs( srcDofIdx ) =
                                     srcDofMemory[volumedofspace::indexing::index( neighborInfo.neighborElementIndices( n ).x(),
                                                                                   neighborInfo.neighborElementIndices( n ).y(),
                                                                                   neighborInfo.neighborFaceType( n ),
                                                                                   srcDofIdx,
                                                                                   numSrcDofs,
                                                                                   level,
                                                                                   srcMemLayout )];
                              }
                              else
                              {
                                 nSrcDofs( srcDofIdx ) =
                                     srcDofMemory[volumedofspace::indexing::index( neighborInfo.neighborElementIndices( n ).x(),
                                                                                   neighborInfo.neighborElementIndices( n ).y(),
                                                                                   neighborInfo.neighborElementIndices( n ).z(),
                                                                                   neighborInfo.neighborCellType( n ),
                                                                                   srcDofIdx,
                                                                                   numSrcDofs,
                                                                                   level,
                                                                                   srcMemLayout )];
                              }
                           }

                           if ( mat == nullptr )
                           {
                              // Matrix-vector multiplication.
                              dstDofs += localMat * nSrcDofs;
                           }
                           else
                           {
                              // TODO: improve this monster
                              std::map< facedof::FaceType, uint_t > invFaceTypeMap;
                              std::map< celldof::CellType, uint_t > invCellTypeMap;

                              for ( uint_t i = 0; i < 2; i++ )
                              {
                                 invFaceTypeMap[facedof::allFaceTypes[i]] = i;
                              }
                              for ( uint_t i = 0; i < 6; i++ )
                              {
                                 invCellTypeMap[celldof::allCellTypes[i]] = i;
                              }

                              uint_t neighborMicroVolType;
                              if ( dim == 2 )
                              {
                                 neighborMicroVolType = invFaceTypeMap[neighborInfo.neighborFaceType( n )];
                              }
                              else
                              {
                                 neighborMicroVolType = invCellTypeMap[neighborInfo.neighborCellType( n )];
                              }
                              addLocalToGlobalMatrix( dim,
                                                      numSrcDofs,
                                                      numDstDofs,
                                                      srcDofMemory,
                                                      dstDofMemory,
                                                      srcMemLayout,
                                                      dstMemLayout,
                                                      neighborInfo.neighborElementIndices( n ),
                                                      elementIdx,
                                                      neighborMicroVolType,
                                                      microVolType,
                                                      level,
                                                      mat,
                                                      localMat );
                           }
                        }
                     }
                  } // End loop over neighboring volumes.
               }    // End if( !onlyVolumeIntegrals() )

               if ( mat == nullptr )
               {
                  // Write DoFs.
                  for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
                  {
                     if ( dim == 2 )
                     {
                        if ( updateType == Replace )
                        {
                           dstDofMemory[volumedofspace::indexing::index(
                               elementIdx.x(), elementIdx.y(), faceType, dstDofIdx, numDstDofs, level, dstMemLayout )] =
                               dstDofs( dstDofIdx );
                        }
                        else if ( updateType == Add )
                        {
                           dstDofMemory[volumedofspace::indexing::index(
                               elementIdx.x(), elementIdx.y(), faceType, dstDofIdx, numDstDofs, level, dstMemLayout )] +=
                               dstDofs( dstDofIdx );
                        }
                        else
                        {
                           WALBERLA_ABORT( "Invalid update type." );
                        }
                     }
                     else
                     {
                        if ( updateType == Replace )
                        {
                           dstDofMemory[volumedofspace::indexing::index( elementIdx.x(),
                                                                         elementIdx.y(),
                                                                         elementIdx.z(),
                                                                         cellType,
                                                                         dstDofIdx,
                                                                         numDstDofs,
                                                                         level,
                                                                         dstMemLayout )] = dstDofs( dstDofIdx );
                        }
                        else if ( updateType == Add )
                        {
                           dstDofMemory[volumedofspace::indexing::index( elementIdx.x(),
                                                                         elementIdx.y(),
                                                                         elementIdx.z(),
                                                                         cellType,
                                                                         dstDofIdx,
                                                                         numDstDofs,
                                                                         level,
                                                                         dstMemLayout )] += dstDofs( dstDofIdx );
                        }
                        else
                        {
                           WALBERLA_ABORT( "Invalid update type." );
                        }
                     }
                  }
               }
            }
         }
      }

      WALBERLA_UNUSED( flag );
   }

   std::shared_ptr< DGForm > form_;
};

} // namespace dg
} // namespace hyteg