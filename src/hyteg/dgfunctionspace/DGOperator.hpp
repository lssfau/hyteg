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
            // WALBERLA_LOG_INFO_ON_ROOT( "Adding value " << localMat( dstDofIdx, srcDofIdx ) << " in ( " << globalRowIdx << ", " << globalColIdx << " )" )
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
         for ( const auto& faceIt : this->getStorage()->getFaces() )
         {
            const auto faceId = faceIt.first;
            const auto face   = *faceIt.second;

            const auto srcPolyDegree = src.polynomialDegree( faceId );
            const auto dstPolyDegree = dst.polynomialDegree( faceId );

            const auto numSrcDofs = src.basis()->numDoFsPerElement( srcPolyDegree );
            const auto numDstDofs = dst.basis()->numDoFsPerElement( dstPolyDegree );

            const auto srcDofMemory = src.volumeDoFFunction()->dofMemory( faceId, level );
            auto       dstDofMemory = dst.volumeDoFFunction()->dofMemory( faceId, level );

            std::map< uint_t, VType* > glMemory;
            for ( const auto& [n, _] : face.getIndirectNeighborFaceIDs() )
            {
               glMemory[n] = src.volumeDoFFunction()->glMemory( faceId, level, n );
            }

            const auto srcMemLayout = src.volumeDoFFunction()->memoryLayout();
            const auto dstMemLayout = dst.volumeDoFFunction()->memoryLayout();

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
               {
                  const auto vertexIndicesVolume = facedof::macroface::getMicroVerticesFromMicroFace( elementIdx, faceType );

                  std::array< Eigen::Matrix< real_t, 2, 1 >, 3 > vertexCoordsVolume;
                  for ( uint_t i = 0; i < 3; i++ )
                  {
                     const auto coord = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndicesVolume[i] );
                     vertexCoordsVolume[i]( 0 ) = coord[0];
                     vertexCoordsVolume[i]( 1 ) = coord[1];
                  }

                  // We only write to the DoFs in the current volume, let's prepare a temporary vector for that.
                  Eigen::Matrix< real_t, Eigen::Dynamic, 1 > dstDofs;
                  dstDofs.resize( numDstDofs, Eigen::NoChange_t::NoChange );
                  dstDofs.setZero();

                  /////////////////////////
                  // Volume contribution //
                  /////////////////////////

                  Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > localMat;
                  form_->integrateVolume(
                      vertexCoordsVolume, *src.basis(), *dst.basis(), srcPolyDegree, dstPolyDegree, localMat );

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

                  {
                     // TODO: all these coord computations can be executed _once_ and then the coordinates can be incremented by h
                     // TODO: blending

                     std::array< Index, 3 >                                          neighborElementIndices;
                     std::array< std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >, 3 > neighborElementVertexCoords;
                     std::array< std::array< Index, 2 >, 3 >                         interfaceVertexIndices;
                     std::array< std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >, 3 > interfaceVertexCoords;
                     std::array< Index, 3 >                                          oppositeVertexIndex;
                     std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >                  oppositeVertexCoords;
                     std::array< Index, 3 >                                          neighborOppositeVertexIndex;
                     std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >                  neighborOppositeVertexCoords;

                     facedof::FaceType nFaceType;

                     if ( faceType == facedof::FaceType::GRAY )
                     {
                        nFaceType = facedof::FaceType::BLUE;

                        neighborElementIndices[0] = Index( elementIdx.x(), elementIdx.y() - 1, 0 );
                        neighborElementIndices[1] = Index( elementIdx.x() - 1, elementIdx.y(), 0 );

                        neighborElementIndices[2] = Index( elementIdx.x(), elementIdx.y(), 0 );

                        // bottom neighbor
                        interfaceVertexIndices[0][0]   = Index( elementIdx.x(), elementIdx.y(), 0 );
                        interfaceVertexIndices[0][1]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
                        oppositeVertexIndex[0]         = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
                        neighborOppositeVertexIndex[0] = Index( elementIdx.x() + 1, elementIdx.y() - 1, 0 );

                        // left neighbor
                        interfaceVertexIndices[1][0]   = Index( elementIdx.x(), elementIdx.y(), 0 );
                        interfaceVertexIndices[1][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
                        oppositeVertexIndex[1]         = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
                        neighborOppositeVertexIndex[1] = Index( elementIdx.x() - 1, elementIdx.y() + 1, 0 );

                        // diagonal neighbor
                        interfaceVertexIndices[2][0]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
                        interfaceVertexIndices[2][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
                        oppositeVertexIndex[2]         = Index( elementIdx.x(), elementIdx.y(), 0 );
                        neighborOppositeVertexIndex[2] = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
                     }
                     else
                     {
                        nFaceType = facedof::FaceType::GRAY;

                        neighborElementIndices[0] = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
                        neighborElementIndices[1] = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
                        neighborElementIndices[2] = Index( elementIdx.x(), elementIdx.y(), 0 );

                        // right neighbor
                        interfaceVertexIndices[0][0]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
                        interfaceVertexIndices[0][1]   = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
                        oppositeVertexIndex[0]         = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
                        neighborOppositeVertexIndex[0] = Index( elementIdx.x() + 2, elementIdx.y(), 0 );

                        // top neighbor
                        interfaceVertexIndices[1][0]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
                        interfaceVertexIndices[1][1]   = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
                        oppositeVertexIndex[1]         = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
                        neighborOppositeVertexIndex[1] = Index( elementIdx.x(), elementIdx.y() + 2, 0 );

                        // diagonal neighbor
                        interfaceVertexIndices[2][0]   = Index( elementIdx.x() + 1, elementIdx.y(), 0 );
                        interfaceVertexIndices[2][1]   = Index( elementIdx.x(), elementIdx.y() + 1, 0 );
                        oppositeVertexIndex[2]         = Index( elementIdx.x() + 1, elementIdx.y() + 1, 0 );
                        neighborOppositeVertexIndex[2] = Index( elementIdx.x(), elementIdx.y(), 0 );
                     }

                     std::array< bool, 3 >    onMacroBoundary;
                     std::array< DoFType, 3 > neighborType;

                     onMacroBoundary[0] = faceType == facedof::FaceType::GRAY && elementIdx.y() == 0;
                     onMacroBoundary[1] = faceType == facedof::FaceType::GRAY && elementIdx.x() == 0;
                     onMacroBoundary[2] =
                         faceType == facedof::FaceType::GRAY &&
                         elementIdx.x() + elementIdx.y() == int( levelinfo::num_microedges_per_edge( level ) - 1 );

                     for ( uint_t n = 0; n < 3; n++ )
                     {
                        neighborType[n] = src.getBoundaryCondition().getBoundaryType(
                            storage_->getEdge( face.neighborEdges()[n] )->getMeshBoundaryFlag() );
                     }

                     // TODO: handle Neumann

                     // Looping over neighbor elements.
                     for ( uint_t n = 0; n < 3; n++ )
                     {
                        const auto vertexIndices =
                            facedof::macroface::getMicroVerticesFromMicroFace( neighborElementIndices[n], nFaceType );

                        for ( uint_t i = 0; i < 3; i++ )
                        {
                           const auto coord = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                           neighborElementVertexCoords[n][i]( 0 ) = coord[0];
                           neighborElementVertexCoords[n][i]( 1 ) = coord[1];
                        }

                        for ( uint_t i = 0; i < 2; i++ )
                        {
                           const auto coord =
                               vertexdof::macroface::coordinateFromIndex( level, face, interfaceVertexIndices[n][i] );
                           interfaceVertexCoords[n][i]( 0 ) = coord[0];
                           interfaceVertexCoords[n][i]( 1 ) = coord[1];
                        }

                        const auto oppCoord = vertexdof::macroface::coordinateFromIndex( level, face, oppositeVertexIndex[n] );
                        oppositeVertexCoords[n]( 0 ) = oppCoord[0];
                        oppositeVertexCoords[n]( 1 ) = oppCoord[1];

                        const auto nOppCoord =
                            vertexdof::macroface::coordinateFromIndex( level, face, neighborOppositeVertexIndex[n] );
                        neighborOppositeVertexCoords[n]( 0 ) = nOppCoord[0];
                        neighborOppositeVertexCoords[n]( 1 ) = nOppCoord[1];

                        // TODO: outsource the normal computation!
                        Eigen::Matrix< real_t, 2, 1 > outerPoint =
                            ( 1 / 3. ) * ( neighborElementVertexCoords[n][0] + neighborElementVertexCoords[n][1] +
                                           neighborElementVertexCoords[n][2] );
                        const auto s = ( outerPoint - interfaceVertexCoords[n][0] )
                                           .template dot( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] ) /
                                       ( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] )
                                           .template dot( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] );
                        const auto proj =
                            interfaceVertexCoords[n][0] + s * ( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] );
                        Eigen::Matrix< real_t, 2, 1 > outwardNormal = ( outerPoint - proj );
                        outwardNormal.normalize();

                        /////////////////////
                        // Domain boundary //
                        /////////////////////

                        if ( onMacroBoundary[n] && neighborType[n] == DirichletBoundary )
                        {
                           ////////////////////////
                           // Dirichlet boundary //
                           ////////////////////////

                           localMat.setZero();
                           form_->integrateFacetDirichletBoundary( vertexCoordsVolume,
                                                                   interfaceVertexCoords[n],
                                                                   oppositeVertexCoords[n],
                                                                   outwardNormal,
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
                        else if ( onMacroBoundary[n] && neighborType[n] == NeumannBoundary )
                        {
                           WALBERLA_ABORT( "Neumann boundary handling not implemented." );
                        }
                        else if ( onMacroBoundary[n] && neighborType[n] == FreeslipBoundary )
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
                           form_->integrateFacetInner( vertexCoordsVolume,
                                                       interfaceVertexCoords[n],
                                                       oppositeVertexCoords[n],
                                                       outwardNormal,
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

                           if ( onMacroBoundary[n] && neighborType[n] == Inner )
                           {
                              ////////////////////////////////////////////////
                              // i) micro-interface on macro-macro-boundary //
                              ////////////////////////////////////////////////

                              // Overwriting some stuff before computing the element matrix.

                              // TODO: some of what follows can be precomputed

                              // The only complicated thing we have to do here is to compute the micro-vertex coordinates
                              // of the micro-vertex that is inside the neighboring macro-volume. Then we must insert the
                              // micro-element coordinates in the _same_ order as this is done locally.

                              // The idea is to perform a "trafo" of the index space to obtain the local logical micro-element
                              // index in the neighboring macro-volume.

                              WALBERLA_ASSERT_GREATER(
                                  face.getIndirectNeighborFaceIDs().count( n ), 0, "Neighborface should exist but doesn't ..." );
                              const auto neighborFace = storage_->getFace( face.getIndirectNeighborFaceIDs().at( n ) );

                              // PID of the opposite macro-vertex.
                              const auto oppositeMacroVertexID =
                                  neighborFace->get_vertex_opposite_to_edge( face.neighborEdges()[n] );

                              Index                   pseudoLocalIndex;
                              std::array< uint_t, 4 > srcBasis;

                              switch ( n )
                              {
                              case 0:
                                 pseudoLocalIndex = Index( elementIdx.x(), 0, 0 );
                                 srcBasis[0]      = neighborFace->vertex_index( face.neighborVertices()[0] );
                                 srcBasis[1]      = neighborFace->vertex_index( face.neighborVertices()[1] );
                                 srcBasis[2]      = neighborFace->vertex_index( oppositeMacroVertexID );
                                 srcBasis[3]      = 3;
                                 break;
                              case 1:
                                 pseudoLocalIndex = Index( elementIdx.y(), 0, 0 );
                                 srcBasis[0]      = neighborFace->vertex_index( face.neighborVertices()[0] );
                                 srcBasis[1]      = neighborFace->vertex_index( face.neighborVertices()[2] );
                                 srcBasis[2]      = neighborFace->vertex_index( oppositeMacroVertexID );
                                 srcBasis[3]      = 3;
                                 break;
                              case 2:
                                 pseudoLocalIndex = Index( elementIdx.y(), 0, 0 );
                                 srcBasis[0]      = neighborFace->vertex_index( face.neighborVertices()[1] );
                                 srcBasis[1]      = neighborFace->vertex_index( face.neighborVertices()[2] );
                                 srcBasis[2]      = neighborFace->vertex_index( oppositeMacroVertexID );
                                 srcBasis[3]      = 3;
                                 break;
                              default:
                                 WALBERLA_ABORT( "Invalid neighbor ID." );
                              }

                              const auto elementIdxNeighborMicro = indexing::basisConversion(
                                  pseudoLocalIndex, srcBasis, { 0, 1, 2, 3 }, levelinfo::num_microedges_per_edge( level ) );

                              // In 2D the neighboring micro-face is always of type GRAY since only those are on the boundary.
                              const auto nFaceTypeLocalToNeighbor = facedof::FaceType::GRAY;

                              const auto neighborElementVertexIndices = facedof::macroface::getMicroVerticesFromMicroFace(
                                  elementIdxNeighborMicro, nFaceTypeLocalToNeighbor );

                              // Eventually we need to know the coordinates of the micro-vertex that is opposite to the interface.
                              // First find out what the local interface ID of the other macro-volume is.
                              uint_t neighborInterfaceID;
                              for ( const auto& [nifID, pid] : neighborFace->getIndirectNeighborFaceIDs() )
                              {
                                 if ( pid == faceId )
                                 {
                                    neighborInterfaceID = nifID;
                                    break;
                                 }
                              }

                              // Then find the idx of the opposing micro-vertex.
                              Index opposingMicroVertexIdx;
                              for ( const auto& nev : neighborElementVertexIndices )
                              {
                                 const auto localInterfaceIDs = isOnCellEdge( nev, level );
                                 if ( !algorithms::contains( localInterfaceIDs, neighborInterfaceID ) )
                                 {
                                    opposingMicroVertexIdx = nev;
                                    break;
                                 }
                              }

                              for ( uint_t i = 0; i < 3; i++ )
                              {
                                 const auto coord = vertexdof::macroface::coordinateFromIndex(
                                     level, *neighborFace, neighborElementVertexIndices[i] );
                                 neighborElementVertexCoords[n][i]( 0 ) = coord[0];
                                 neighborElementVertexCoords[n][i]( 1 ) = coord[1];
                              }

                              const auto coord =
                                  vertexdof::macroface::coordinateFromIndex( level, *neighborFace, opposingMicroVertexIdx );

                              neighborOppositeVertexCoords[n]( 0 ) = coord[0];
                              neighborOppositeVertexCoords[n]( 1 ) = coord[1];

                              // We need to recompute the normal (?)...
                              // TODO: outsource the normal computation!
                              Eigen::Matrix< real_t, 2, 1 > outerPoint =
                                  ( 1 / 3. ) * ( neighborElementVertexCoords[n][0] + neighborElementVertexCoords[n][1] +
                                                 neighborElementVertexCoords[n][2] );
                              const auto s = ( outerPoint - interfaceVertexCoords[n][0] )
                                                 .template dot( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] ) /
                                             ( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] )
                                                 .template dot( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] );
                              const auto proj =
                                  interfaceVertexCoords[n][0] + s * ( interfaceVertexCoords[n][1] - interfaceVertexCoords[n][0] );
                              Eigen::Matrix< real_t, 2, 1 > outwardNormal = ( outerPoint - proj );
                              outwardNormal.normalize();

                              localMat.setZero();
                              form_->integrateFacetCoupling( vertexCoordsVolume,
                                                             neighborElementVertexCoords[n],
                                                             interfaceVertexCoords[n],
                                                             oppositeVertexCoords[n],
                                                             neighborOppositeVertexCoords[n],
                                                             outwardNormal,
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
                                 nSrcDoFArrIndices[srcDofIdx] =
                                     volumedofspace::indexing::indexGhostLayer( n,
                                                                                neighborElementIndices[n].x(),
                                                                                neighborElementIndices[n].y(),
                                                                                nFaceType,
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
                                       const auto globalRowIdx = dstDofMemory[volumedofspace::indexing::index( elementIdx.x(),
                                                                                                               elementIdx.y(),
                                                                                                               faceType,
                                                                                                               dstDofIdx,
                                                                                                               numDstDofs,
                                                                                                               level,
                                                                                                               dstMemLayout )];
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
                              form_->integrateFacetCoupling( vertexCoordsVolume,
                                                             neighborElementVertexCoords[n],
                                                             interfaceVertexCoords[n],
                                                             oppositeVertexCoords[n],
                                                             neighborOppositeVertexCoords[n],
                                                             outwardNormal,
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
                                     srcDofMemory[volumedofspace::indexing::index( neighborElementIndices[n].x(),
                                                                                   neighborElementIndices[n].y(),
                                                                                   nFaceType,
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
                                                         neighborElementIndices[n],
                                                         elementIdx,
                                                         nFaceType,
                                                         faceType,
                                                         level,
                                                         mat,
                                                         localMat );
                              }
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