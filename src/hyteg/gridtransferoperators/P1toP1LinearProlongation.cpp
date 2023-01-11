/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_2D_macroface_P1_push_additive.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_3D_macrocell_P1_push_additive.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

namespace hyteg {

void P1toP1LinearProlongation::prolongate2D( const P1Function< real_t >& function,
                                             const uint_t&               sourceLevel,
                                             const DoFType&              flag ) const
{
   const uint_t destinationLevel = sourceLevel + 1;

   function.communicate< Vertex, Edge >( sourceLevel );
   function.communicate< Edge, Face >( sourceLevel );
   function.communicate< Face, Edge >( sourceLevel );
   function.communicate< Edge, Vertex >( sourceLevel );

   for ( const auto& it : function.getStorage()->getVertices() )
   {
      const Vertex& vertex = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         const auto srcData = vertex.getData( function.getVertexDataID() )->getPointer( sourceLevel );
         auto       dstData = vertex.getData( function.getVertexDataID() )->getPointer( destinationLevel );
         prolongateMacroVertex2D( srcData, dstData, sourceLevel );
      }
   }

   for ( const auto& it : function.getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         const auto srcData = edge.getData( function.getEdgeDataID() )->getPointer( sourceLevel );
         auto       dstData = edge.getData( function.getEdgeDataID() )->getPointer( destinationLevel );
         prolongateMacroEdge2D( srcData, dstData, sourceLevel );
      }
   }

   for ( const auto& it : function.getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      if ( testFlag( function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         const auto srcData = face.getData( function.getFaceDataID() )->getPointer( sourceLevel );
         auto       dstData = face.getData( function.getFaceDataID() )->getPointer( destinationLevel );
         prolongateMacroFace2D( srcData, dstData, sourceLevel );
      }
   }
}

static real_t calculateInverseFactorToScaleNeighborhoodContribution( const std::array< real_t, 4 >& invNumNeighborsOfVertex,
                                                                     const std::array< real_t, 6 >& invNumNeighborsOfEdge,
                                                                     const std::array< real_t, 4 >  invNumNeighborsOfFace,
                                                                     const indexing::Index&         index,
                                                                     const uint_t&                  level )
{
   const auto stencilLeafOnCellVertices = vertexdof::macrocell::isOnCellVertex( index, level );
   const auto stencilLeafOnCellEdges    = vertexdof::macrocell::isOnCellEdge( index, level );
   const auto stencilLeafOnCellFaces    = vertexdof::macrocell::isOnCellFace( index, level );

   real_t invFactorDueToNeighborhood = real_c( 1 );

   if ( stencilLeafOnCellVertices.size() > 0 )
   {
      WALBERLA_ASSERT_EQUAL( stencilLeafOnCellVertices.size(), 1 );
      const auto localVertexID   = *stencilLeafOnCellVertices.begin();
      invFactorDueToNeighborhood = invNumNeighborsOfVertex[localVertexID];
   }
   else if ( stencilLeafOnCellEdges.size() > 0 )
   {
      WALBERLA_ASSERT_EQUAL( stencilLeafOnCellEdges.size(), 1 );
      const auto localEdgeID     = *stencilLeafOnCellEdges.begin();
      invFactorDueToNeighborhood = invNumNeighborsOfEdge[localEdgeID];
   }
   else if ( stencilLeafOnCellFaces.size() > 0 )
   {
      WALBERLA_ASSERT_EQUAL( stencilLeafOnCellFaces.size(), 1 );
      const auto localFaceID     = *stencilLeafOnCellFaces.begin();
      invFactorDueToNeighborhood = invNumNeighborsOfFace[localFaceID];
   }

   return invFactorDueToNeighborhood;
}

void P1toP1LinearProlongation::prolongate2DAdditively( const P1Function< real_t >& function,
                                                       const uint_t&               sourceLevel,
                                                       const DoFType&              flag,
                                                       const UpdateType&           updateType ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
   const DoFType excludeFlag = ( flag ^ All );

   const uint_t destinationLevel = sourceLevel + 1;

   function.communicate< Vertex, Edge >( sourceLevel );
   function.communicate< Edge, Face >( sourceLevel );

   auto storage = function.getStorage();

   for ( const auto& faceIt : function.getStorage()->getFaces() )
   {
      const auto face    = faceIt.second;
      const auto srcData = face->getData( function.getFaceDataID() )->getPointer( sourceLevel );
      auto       dstData = face->getData( function.getFaceDataID() )->getPointer( destinationLevel );

      if ( updateType == Replace )
      {
         // Clear complete face (incl. ghost-layer) on dst level
         for ( const auto& dstIdx : vertexdof::macroface::Iterator( destinationLevel, 0 ) )
         {
            const auto arrayIdxDst = vertexdof::macroface::index( destinationLevel, dstIdx.x(), dstIdx.y() );
            dstData[arrayIdxDst]   = real_c( 0 );
         }
      }
      else if ( updateType == Add )
      {
         // we only set the ghost layers to zero, but not the inner unknowns
         for ( const auto& dstIdx : vertexdof::macroface::Iterator( destinationLevel, 0 ) )
         {
            if ( vertexdof::macroface::isVertexOnBoundary( destinationLevel, dstIdx ) )
            {
               const auto arrayIdxDst = vertexdof::macroface::index( destinationLevel, dstIdx.x(), dstIdx.y() );
               dstData[arrayIdxDst]   = real_c( 0 );
            }
         }
      }
      else
      {
         WALBERLA_ABORT( "Invalid update type in prolongation." );
      }

      const auto numNeighborFacesEdge0 =
          static_cast< real_t >( storage->getEdge( face->neighborEdges().at( 0 ) )->getNumNeighborFaces() );
      const auto numNeighborFacesEdge1 =
          static_cast< real_t >( storage->getEdge( face->neighborEdges().at( 1 ) )->getNumNeighborFaces() );
      const auto numNeighborFacesEdge2 =
          static_cast< real_t >( storage->getEdge( face->neighborEdges().at( 2 ) )->getNumNeighborFaces() );
      const auto numNeighborFacesVertex0 =
          static_cast< real_t >( storage->getVertex( face->neighborVertices().at( 0 ) )->getNumNeighborFaces() );
      const auto numNeighborFacesVertex1 =
          static_cast< real_t >( storage->getVertex( face->neighborVertices().at( 1 ) )->getNumNeighborFaces() );
      const auto numNeighborFacesVertex2 =
          static_cast< real_t >( storage->getVertex( face->neighborVertices().at( 2 ) )->getNumNeighborFaces() );

      vertexdof::macroface::generated::prolongate_2D_macroface_P1_push_additive( srcData,
                                                                                 dstData,
                                                                                 static_cast< int32_t >( sourceLevel ),
                                                                                 numNeighborFacesEdge0,
                                                                                 numNeighborFacesEdge1,
                                                                                 numNeighborFacesEdge2,
                                                                                 numNeighborFacesVertex0,
                                                                                 numNeighborFacesVertex1,
                                                                                 numNeighborFacesVertex2 );
   }

   function.communicateAdditively< Face, Edge >( destinationLevel, excludeFlag, *storage, updateType == Replace );
   function.communicateAdditively< Face, Vertex >( destinationLevel, excludeFlag, *storage, updateType == Replace );
}

void P1toP1LinearProlongation::prolongate3DAdditively( const P1Function< real_t >& function,
                                                       const uint_t&               sourceLevel,
                                                       const DoFType&              flag,
                                                       const UpdateType&           updateType ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
   const DoFType excludeFlag = ( flag ^ All );

   const uint_t destinationLevel = sourceLevel + 1;

   function.communicate< Vertex, Edge >( sourceLevel );
   function.communicate< Edge, Face >( sourceLevel );
   function.communicate< Face, Cell >( sourceLevel );

   for ( const auto& cellIt : function.getStorage()->getCells() )
   {
      const auto cell    = cellIt.second;
      const auto srcData = cell->getData( function.getCellDataID() )->getPointer( sourceLevel );
      auto       dstData = cell->getData( function.getCellDataID() )->getPointer( destinationLevel );

      if ( updateType == Replace )
      {
         // Clear complete face (incl. ghost-layer) on dst level
         for ( const auto& dstIdx : vertexdof::macrocell::Iterator( destinationLevel, 0 ) )
         {
            const auto arrayIdxDst = vertexdof::macrocell::index( destinationLevel, dstIdx.x(), dstIdx.y(), dstIdx.z() );
            dstData[arrayIdxDst]   = real_c( 0 );
         }
      }
      else if ( updateType == Add )
      {
         // we only set the ghost layers to zero, but not the inner unknowns
         for ( const auto& dstIdx : vertexdof::macrocell::Iterator( destinationLevel, 0 ) )
         {
            if ( !vertexdof::macrocell::isOnCellFace( dstIdx, destinationLevel ).empty() )
            {
               const auto arrayIdxDst = vertexdof::macrocell::index( destinationLevel, dstIdx.x(), dstIdx.y(), dstIdx.z() );
               dstData[arrayIdxDst]   = real_c( 0 );
            }
         }
      }
      else
      {
         WALBERLA_ABORT( "Invalid update type in prolongation." );
      }

      if ( globalDefines::useGeneratedKernels )
      {
#ifdef WALBERLA_USE_GENERATED_KERNELS
         auto storage = function.getStorage();

         const double numNeighborCellsFace0 =
             static_cast< double >( storage->getFace( cell->neighborFaces().at( 0 ) )->getNumNeighborCells() );
         const double numNeighborCellsFace1 =
             static_cast< double >( storage->getFace( cell->neighborFaces().at( 1 ) )->getNumNeighborCells() );
         const double numNeighborCellsFace2 =
             static_cast< double >( storage->getFace( cell->neighborFaces().at( 2 ) )->getNumNeighborCells() );
         const double numNeighborCellsFace3 =
             static_cast< double >( storage->getFace( cell->neighborFaces().at( 3 ) )->getNumNeighborCells() );

         const double numNeighborCellsEdge0 =
             static_cast< double >( storage->getEdge( cell->neighborEdges().at( 0 ) )->getNumNeighborCells() );
         const double numNeighborCellsEdge1 =
             static_cast< double >( storage->getEdge( cell->neighborEdges().at( 1 ) )->getNumNeighborCells() );
         const double numNeighborCellsEdge2 =
             static_cast< double >( storage->getEdge( cell->neighborEdges().at( 2 ) )->getNumNeighborCells() );
         const double numNeighborCellsEdge3 =
             static_cast< double >( storage->getEdge( cell->neighborEdges().at( 3 ) )->getNumNeighborCells() );
         const double numNeighborCellsEdge4 =
             static_cast< double >( storage->getEdge( cell->neighborEdges().at( 4 ) )->getNumNeighborCells() );
         const double numNeighborCellsEdge5 =
             static_cast< double >( storage->getEdge( cell->neighborEdges().at( 5 ) )->getNumNeighborCells() );

         const double numNeighborCellsVertex0 =
             static_cast< double >( storage->getVertex( cell->neighborVertices().at( 0 ) )->getNumNeighborCells() );
         const double numNeighborCellsVertex1 =
             static_cast< double >( storage->getVertex( cell->neighborVertices().at( 1 ) )->getNumNeighborCells() );
         const double numNeighborCellsVertex2 =
             static_cast< double >( storage->getVertex( cell->neighborVertices().at( 2 ) )->getNumNeighborCells() );
         const double numNeighborCellsVertex3 =
             static_cast< double >( storage->getVertex( cell->neighborVertices().at( 3 ) )->getNumNeighborCells() );

         vertexdof::macrocell::generated::prolongate_3D_macrocell_P1_push_additive( srcData,
                                                                                    dstData,
                                                                                    static_cast< int32_t >( sourceLevel ),
                                                                                    numNeighborCellsEdge0,
                                                                                    numNeighborCellsEdge1,
                                                                                    numNeighborCellsEdge2,
                                                                                    numNeighborCellsEdge3,
                                                                                    numNeighborCellsEdge4,
                                                                                    numNeighborCellsEdge5,
                                                                                    numNeighborCellsFace0,
                                                                                    numNeighborCellsFace1,
                                                                                    numNeighborCellsFace2,
                                                                                    numNeighborCellsFace3,
                                                                                    numNeighborCellsVertex0,
                                                                                    numNeighborCellsVertex1,
                                                                                    numNeighborCellsVertex2,
                                                                                    numNeighborCellsVertex3 );
#endif
      }
      else
      {
         // Calculate inverse number of neighboring cells for each neighboring macro-primitive.
         std::array< real_t, 4 > invNumNeighborsOfVertex;
         std::array< real_t, 6 > invNumNeighborsOfEdge;
         std::array< real_t, 4 > invNumNeighborsOfFace;

         invNumNeighborsOfVertex.fill( real_c( 0 ) );
         invNumNeighborsOfEdge.fill( real_c( 0 ) );
         invNumNeighborsOfFace.fill( real_c( 0 ) );

         for ( const auto& neighborVertexID : cell->neighborVertices() )
         {
            invNumNeighborsOfVertex[cell->getLocalVertexID( neighborVertexID )] =
                real_c( 1 ) / real_c( function.getStorage()->getVertex( neighborVertexID )->getNumNeighborCells() );
         }
         for ( const auto& neighborEdgeID : cell->neighborEdges() )
         {
            invNumNeighborsOfEdge[cell->getLocalEdgeID( neighborEdgeID )] =
                real_c( 1 ) / real_c( function.getStorage()->getEdge( neighborEdgeID )->getNumNeighborCells() );
         }
         for ( const auto& neighborFaceID : cell->neighborFaces() )
         {
            invNumNeighborsOfFace[cell->getLocalFaceID( neighborFaceID )] =
                real_c( 1 ) / real_c( function.getStorage()->getFace( neighborFaceID )->getNumNeighborCells() );
         }

         for ( const auto& srcIdx : vertexdof::macrocell::Iterator( sourceLevel ) )
         {
            const auto arrayIdxSrc = vertexdof::macrocell::index( sourceLevel, srcIdx.x(), srcIdx.y(), srcIdx.z() );
            const auto dstIdx      = srcIdx * 2;

            const auto onCellVertices = vertexdof::macrocell::isOnCellVertex( srcIdx, sourceLevel );
            const auto onCellEdges    = vertexdof::macrocell::isOnCellEdge( srcIdx, sourceLevel );
            const auto onCellFaces    = vertexdof::macrocell::isOnCellFace( srcIdx, sourceLevel );

            // update center
            const auto invFactorToScaleContributionCenter = calculateInverseFactorToScaleNeighborhoodContribution(
                invNumNeighborsOfVertex, invNumNeighborsOfEdge, invNumNeighborsOfFace, dstIdx, destinationLevel );

            const auto arrayIdxDstCenter = vertexdof::macrocell::index( destinationLevel, dstIdx.x(), dstIdx.y(), dstIdx.z() );
            dstData[arrayIdxDstCenter] += invFactorToScaleContributionCenter * srcData[arrayIdxSrc];

            // update new points depending on location in macro-cell
            if ( onCellVertices.size() > 0 )
            {
               WALBERLA_ASSERT_EQUAL( onCellVertices.size(), 1 );
               const auto localVertexID = *onCellVertices.begin();

               for ( const auto& dir : vertexdof::macrocell::neighborsOnVertexWithoutCenter[localVertexID] )
               {
                  const auto increment                    = vertexdof::logicalIndexOffsetFromVertex( dir );
                  const auto dirIdxDst                    = dstIdx + increment;
                  const auto invFactorToScaleContribution = calculateInverseFactorToScaleNeighborhoodContribution(
                      invNumNeighborsOfVertex, invNumNeighborsOfEdge, invNumNeighborsOfFace, dirIdxDst, destinationLevel );

                  const auto arrayIdxDst =
                      vertexdof::macrocell::index( destinationLevel, dirIdxDst.x(), dirIdxDst.y(), dirIdxDst.z() );
                  dstData[arrayIdxDst] += real_c( 0.5 ) * invFactorToScaleContribution * srcData[arrayIdxSrc];
               }
            }
            else if ( onCellEdges.size() > 0 )
            {
               WALBERLA_ASSERT_EQUAL( onCellEdges.size(), 1 );
               const auto localEdgeID = *onCellEdges.begin();

               for ( const auto& dir : vertexdof::macrocell::neighborsOnEdgeWithoutCenter[localEdgeID] )
               {
                  const auto increment                    = vertexdof::logicalIndexOffsetFromVertex( dir );
                  const auto dirIdxDst                    = dstIdx + increment;
                  const auto invFactorToScaleContribution = calculateInverseFactorToScaleNeighborhoodContribution(
                      invNumNeighborsOfVertex, invNumNeighborsOfEdge, invNumNeighborsOfFace, dirIdxDst, destinationLevel );
                  const auto arrayIdxDst =
                      vertexdof::macrocell::index( destinationLevel, dirIdxDst.x(), dirIdxDst.y(), dirIdxDst.z() );
                  dstData[arrayIdxDst] += real_c( 0.5 ) * invFactorToScaleContribution * srcData[arrayIdxSrc];
               }
            }
            else if ( onCellFaces.size() > 0 )
            {
               WALBERLA_ASSERT_EQUAL( onCellFaces.size(), 1 );
               const auto localFaceID = *onCellFaces.begin();

               for ( const auto& dir : vertexdof::macrocell::neighborsOnFaceWithoutCenter[localFaceID] )
               {
                  const auto increment                    = vertexdof::logicalIndexOffsetFromVertex( dir );
                  const auto dirIdxDst                    = dstIdx + increment;
                  const auto invFactorToScaleContribution = calculateInverseFactorToScaleNeighborhoodContribution(
                      invNumNeighborsOfVertex, invNumNeighborsOfEdge, invNumNeighborsOfFace, dirIdxDst, destinationLevel );
                  const auto arrayIdxDst =
                      vertexdof::macrocell::index( destinationLevel, dirIdxDst.x(), dirIdxDst.y(), dirIdxDst.z() );
                  dstData[arrayIdxDst] += real_c( 0.5 ) * invFactorToScaleContribution * srcData[arrayIdxSrc];
               }
            }
            else
            {
               for ( const auto& dir : vertexdof::macrocell::neighborsWithoutCenter )
               {
                  const auto arrayIdxDst =
                      vertexdof::macrocell::indexFromVertex( destinationLevel, dstIdx.x(), dstIdx.y(), dstIdx.z(), dir );
                  dstData[arrayIdxDst] += real_c( 0.5 ) * srcData[arrayIdxSrc];
               }
            }
         }
      }
   }

   function.communicateAdditively< Cell, Vertex >( destinationLevel, excludeFlag, *function.getStorage(), updateType == Replace );
   function.communicateAdditively< Cell, Edge >( destinationLevel, excludeFlag, *function.getStorage(), updateType == Replace );
   function.communicateAdditively< Cell, Face >( destinationLevel, excludeFlag, *function.getStorage(), updateType == Replace );
}

void P1toP1LinearProlongation::prolongateMacroVertex2D( const real_t* src, real_t* dst, const uint_t& ) const
{
   dst[0] = src[0];
}

void P1toP1LinearProlongation::prolongateMacroEdge2D( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   uint_t rowsize_c = levelinfo::num_microvertices_per_edge( sourceLevel );
   idx_t i_c;

   for ( i_c = 1; i_c < idx_t( rowsize_c ) - 1; ++i_c )
   {
      dst[vertexdof::macroedge::indexFromVertex( sourceLevel + 1, 2 * i_c, stencilDirection::VERTEX_C )] =
          src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c, stencilDirection::VERTEX_C )];

      dst[vertexdof::macroedge::indexFromVertex( sourceLevel + 1, 2 * i_c - 1, stencilDirection::VERTEX_C )] =
          real_c( 0.5 ) * ( src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c - 1, stencilDirection::VERTEX_C )] +
                            src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c, stencilDirection::VERTEX_C )] );
   }

   dst[vertexdof::macroedge::indexFromVertex( sourceLevel + 1, 2 * i_c - 1, stencilDirection::VERTEX_C )] =
       real_c( 0.5 ) * ( src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c - 1, stencilDirection::VERTEX_C )] +
                         src[vertexdof::macroedge::indexFromVertex( sourceLevel, i_c, stencilDirection::VERTEX_C )] );
}

void P1toP1LinearProlongation::prolongateMacroFace2D( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   typedef stencilDirection SD;
   using namespace vertexdof::macroface;

   uint_t N_c   = levelinfo::num_microvertices_per_edge( sourceLevel );
   uint_t N_c_i = N_c;

   idx_t j;

   for ( idx_t i = 1; i < idx_t( N_c ) - 1; ++i )
   {
      for ( j = 1; j < idx_t( N_c_i ) - 2; ++j )
      {
         dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j, SD::VERTEX_C )] =
             src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )];
         dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] =
             real_c( 0.5 ) * ( src[indexFromVertex( sourceLevel, i - 1, j, SD::VERTEX_C )] +
                               src[indexFromVertex( sourceLevel, i, j - 1, SD::VERTEX_C )] );
         dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] =
             real_c( 0.5 ) * ( src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )] +
                               src[indexFromVertex( sourceLevel, i - 1, j, SD::VERTEX_C )] );
         dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] =
             real_c( 0.5 ) * ( src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )] +
                               src[indexFromVertex( sourceLevel, i, j - 1, SD::VERTEX_C )] );
      }

      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] =
          real_c( 0.5 ) * ( src[indexFromVertex( sourceLevel, i - 1, j, SD::VERTEX_C )] +
                            src[indexFromVertex( sourceLevel, i, j - 1, SD::VERTEX_C )] );
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] =
          real_c( 0.5 ) * ( src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )] +
                            src[indexFromVertex( sourceLevel, i - 1, j, SD::VERTEX_C )] );
      dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] =
          real_c( 0.5 ) * ( src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )] +
                            src[indexFromVertex( sourceLevel, i, j - 1, SD::VERTEX_C )] );

      --N_c_i;
   }
}

} // namespace hyteg
