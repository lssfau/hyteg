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

using facedof::FaceType;
using indexing::Index;
using volumedofspace::indexing::VolumeDoFMemoryLayout;
using walberla::int_c;
using walberla::real_t;

class DGOperator : public Operator< DGFunction< real_t >, DGFunction< real_t > >,
                   public WeightedJacobiSmoothable< DGFunction< real_t > >,
                   public OperatorWithInverseDiagonal< DGFunction< real_t > >
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

   void smooth_jac( const DGFunction< real_t >& dst,
                    const DGFunction< real_t >& rhs,
                    const DGFunction< real_t >& tmp,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override;

   [[nodiscard]] std::shared_ptr< DGFunction< real_t > > getInverseDiagonalValues() const override;

 private:
   /// Just a small helper method that writes the local matrix into the global sparse system.
   template < typename VType >
   void addLocalToGlobalMatrix( int                                                            numSrcDofs,
                                int                                                            numDstDofs,
                                VType*                                                         srcDofMemory,
                                VType*                                                         dstDofMemory,
                                VolumeDoFMemoryLayout                                          srcMemLayout,
                                VolumeDoFMemoryLayout                                          dstMemLayout,
                                Index                                                          srcElementIdx,
                                Index                                                          dstElementIdx,
                                FaceType                                                       srcFaceType,
                                FaceType                                                       dstFaceType,
                                uint_t                                                         level,
                                std::shared_ptr< SparseMatrixProxy >                           mat,
                                const Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const
   {
      // Sparse assembly.
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
      using indexing::Index;

      WALBERLA_CHECK( updateType == Replace );

      src.communicate( level );

      if ( this->getStorage()->hasGlobalCells() )
      {
         WALBERLA_ABORT( "DG apply not implemented in 3D." );
      }
      else
      {
         const int dim = 2;

         for ( const auto& faceIt : this->getStorage()->getFaces() )
         {
            const auto  faceId = faceIt.first;
            const auto& face   = *faceIt.second;

            const auto srcPolyDegree = src.polynomialDegree( faceId );
            const auto dstPolyDegree = dst.polynomialDegree( faceId );

            const auto numSrcDofs = src.basis()->numDoFsPerElement( srcPolyDegree );
            const auto numDstDofs = dst.basis()->numDoFsPerElement( dstPolyDegree );

            const auto srcDofMemory = src.volumeDoFFunction()->dofMemory( faceId, level );
            auto       dstDofMemory = dst.volumeDoFFunction()->dofMemory( faceId, level );

            std::map< uint_t, VType* > glMemory;
            for ( const auto& [n, _] : face.getIndirectNeighborFaceIDsOverEdges() )
            {
               glMemory[n] = src.volumeDoFFunction()->glMemory( faceId, level, n );
            }

            const auto srcMemLayout = src.volumeDoFFunction()->memoryLayout();
            const auto dstMemLayout = dst.volumeDoFFunction()->memoryLayout();

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
               {
                  // TODO: all these coord computations can be executed _once_ and then the coordinates can be incremented by h
                  // TODO: blending

                  // This object does the heavy lifting of computing all required coordinates and normals.
                  volumedofspace::indexing::ElementNeighborInfo neighborInfo(
                      elementIdx, faceType, level, src.getBoundaryCondition(), faceId, storage_ );

                  // We only write to the DoFs in the current volume, let's prepare a temporary vector for that.
                  Eigen::Matrix< real_t, Eigen::Dynamic, 1 > dstDofs;
                  dstDofs.resize( numDstDofs, Eigen::NoChange_t::NoChange );
                  dstDofs.setZero();

                  /////////////////////////
                  // Volume contribution //
                  /////////////////////////

                  Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > localMat;
                  form_->integrateVolume( dim,
                                          neighborInfo.elementVertexCoords(),
                                          *src.basis(),
                                          *dst.basis(),
                                          srcPolyDegree,
                                          dstPolyDegree,
                                          localMat );

                  // Volume DoFs are source.
                  Eigen::Matrix< real_t, Eigen::Dynamic, 1 > srcDofs;
                  srcDofs.resize( numSrcDofs, Eigen::NoChange_t::NoChange );

                  for ( uint_t srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
                  {
                     srcDofs( srcDofIdx ) = srcDofMemory[volumedofspace::indexing::index(
                         elementIdx.x(), elementIdx.y(), faceType, srcDofIdx, numSrcDofs, level, srcMemLayout )];
                  }

                  if ( mat == nullptr )
                  {
                     // Matrix-vector multiplication.
                     dstDofs += localMat * srcDofs;
                  }
                  else
                  {
                     // Sparse assembly.
                     addLocalToGlobalMatrix( numSrcDofs,
                                             numDstDofs,
                                             srcDofMemory,
                                             dstDofMemory,
                                             srcMemLayout,
                                             dstMemLayout,
                                             elementIdx,
                                             elementIdx,
                                             faceType,
                                             faceType,
                                             level,
                                             mat,
                                             localMat );
                  }

                  /////////////////////////////
                  // Interface contributions //
                  /////////////////////////////

                  for ( uint_t n = 0; n < 3; n++ )
                  {
                     /////////////////////
                     // Domain boundary //
                     /////////////////////

                     if ( neighborInfo.onMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == DirichletBoundary )
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
                           addLocalToGlobalMatrix( numSrcDofs,
                                                   numDstDofs,
                                                   srcDofMemory,
                                                   dstDofMemory,
                                                   srcMemLayout,
                                                   dstMemLayout,
                                                   elementIdx,
                                                   elementIdx,
                                                   faceType,
                                                   faceType,
                                                   level,
                                                   mat,
                                                   localMat );
                        }
                     }
                     else if ( neighborInfo.onMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == NeumannBoundary )
                     {
                        WALBERLA_ABORT( "Neumann boundary handling not implemented." );
                     }
                     else if ( neighborInfo.onMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == FreeslipBoundary )
                     {
                        WALBERLA_ABORT( "Free-slip boundary handling not implemented." );
                     }

                     //////////////////
                     // Inner domain //
                     //////////////////

                     else
                     {
                        ///////////////////////////////////
                        // a) inner element contribution //
                        ///////////////////////////////////

                        localMat.setZero();
                        form_->integrateFacetInner( dim,
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
                           addLocalToGlobalMatrix( numSrcDofs,
                                                   numDstDofs,
                                                   srcDofMemory,
                                                   dstDofMemory,
                                                   srcMemLayout,
                                                   dstMemLayout,
                                                   elementIdx,
                                                   elementIdx,
                                                   faceType,
                                                   faceType,
                                                   level,
                                                   mat,
                                                   localMat );
                        }

                        ////////////////////////////////////////
                        // b) coupling to neighboring element //
                        ////////////////////////////////////////

                        if ( neighborInfo.onMacroBoundary( n ) && neighborInfo.neighborBoundaryType( n ) == Inner )
                        {
                           ////////////////////////////////////////////////
                           // i) micro-interface on macro-macro-boundary //
                           ////////////////////////////////////////////////

                           // The neighboring micro-element coords have to be computed since they are now different as for an
                           // element on the same macro-volume.
                           std::vector< Eigen::Matrix< real_t, 3, 1 > > neighborElementVertexCoords;
                           Eigen::Matrix< real_t, 3, 1 >                neighborOppositeVertexCoords;

                           neighborInfo.macroBoundaryNeighborElementVertexCoords(
                               n, neighborElementVertexCoords, neighborOppositeVertexCoords );

                           localMat.setZero();
                           form_->integrateFacetCoupling( dim,
                                                          neighborInfo.elementVertexCoords(),
                                                          neighborElementVertexCoords,
                                                          neighborInfo.interfaceVertexCoords( n ),
                                                          neighborInfo.oppositeVertexCoords( n ),
                                                          neighborOppositeVertexCoords,
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
                              // This access might seem a little unintuitive, but it does another bit of lifting.
                              // The ghost-layer data can be accessed with this indexing function by "extending" the macro-volume
                              // structure. One of the indices may now be -1 for example.
                              nSrcDoFArrIndices[srcDofIdx] =
                                  volumedofspace::indexing::indexGhostLayer( n,
                                                                             neighborInfo.neighborElementIndices( n ).x(),
                                                                             neighborInfo.neighborElementIndices( n ).y(),
                                                                             facedof::FaceType::BLUE,
                                                                             srcDofIdx,
                                                                             numSrcDofs,
                                                                             level,
                                                                             srcMemLayout );
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
                                    const auto globalRowIdx = dstDofMemory[volumedofspace::indexing::index(
                                        elementIdx.x(), elementIdx.y(), faceType, dstDofIdx, numDstDofs, level, dstMemLayout )];
                                    const auto globalColIdx = glMemory[n][nSrcDoFArrIndices[srcDofIdx]];
                                    mat->addValue( globalRowIdx, globalColIdx, localMat( dstDofIdx, srcDofIdx ) );
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
                              nSrcDofs( srcDofIdx ) =
                                  srcDofMemory[volumedofspace::indexing::index( neighborInfo.neighborElementIndices( n ).x(),
                                                                                neighborInfo.neighborElementIndices( n ).y(),
                                                                                neighborInfo.neighborFaceType( n ),
                                                                                srcDofIdx,
                                                                                numSrcDofs,
                                                                                level,
                                                                                srcMemLayout )];
                           }

                           if ( mat == nullptr )
                           {
                              // Matrix-vector multiplication.
                              dstDofs += localMat * nSrcDofs;
                           }
                           else
                           {
                              // Sparse assembly.
                              addLocalToGlobalMatrix( numSrcDofs,
                                                      numDstDofs,
                                                      srcDofMemory,
                                                      dstDofMemory,
                                                      srcMemLayout,
                                                      dstMemLayout,
                                                      neighborInfo.neighborElementIndices( n ),
                                                      elementIdx,
                                                      neighborInfo.neighborFaceType( n ),
                                                      faceType,
                                                      level,
                                                      mat,
                                                      localMat );
                           }
                        }
                     }
                  }

                  if ( mat == nullptr )
                  {
                     // Write DoFs.
                     for ( uint_t dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
                     {
                        dstDofMemory[volumedofspace::indexing::index(
                            elementIdx.x(), elementIdx.y(), faceType, dstDofIdx, numDstDofs, level, dstMemLayout )] =
                            dstDofs( dstDofIdx );
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