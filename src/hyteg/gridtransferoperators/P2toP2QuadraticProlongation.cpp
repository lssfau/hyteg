/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"

#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2Multigrid.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_2D_macroface_P2_push_from_vertexdofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_2D_macroface_P2_push_from_edgedofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_2D_macroface_P2_push_from_edgedofs_level_0_to_1.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_3D_macrocell_P2_push_from_vertexdofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_3D_macrocell_P2_push_from_edgedofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/prolongate_3D_macrocell_P2_push_from_edgedofs_level_0_to_1.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

void P2toP2QuadraticProlongation::prolongate( const P2Function< walberla::real_t >& function,
                                              const walberla::uint_t&               sourceLevel,
                                              const DoFType&                        flag ) const
{
   if ( function.getStorage()->hasGlobalCells() )
   {
      prolongateAdditively3D( function, sourceLevel, flag, Replace );
   }
   else
   {
      prolongateAdditively( function, sourceLevel, flag, Replace );
   }
}

void P2toP2QuadraticProlongation::prolongateAndAdd( const P2Function< walberla::real_t >& function,
                                                    const walberla::uint_t&               sourceLevel,
                                                    const DoFType&                        flag ) const
{
   if ( function.getStorage()->hasGlobalCells() )
   {
      prolongateAdditively3D( function, sourceLevel, flag, Add );
   }
   else
   {
      prolongateAdditively( function, sourceLevel, flag, Add );
   }
}

void P2toP2QuadraticProlongation::prolongateAdditively( const P2Function< real_t >& function,
                                                        const uint_t&               sourceLevel,
                                                        const DoFType&              flag,
                                                        const UpdateType &          updateType ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
   const DoFType excludeFlag = ( flag ^ All );

   const auto storage = function.getStorage();

   const uint_t fineLevel   = sourceLevel + 1;
   const uint_t coarseLevel = sourceLevel;

   function.communicate< Vertex, Edge >( coarseLevel );
   function.communicate< Edge, Face >( coarseLevel );

   for ( const auto& faceIt : function.getStorage()->getFaces() )
   {
      const auto face = faceIt.second;

      auto vertexFineData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( fineLevel );
      auto edgeFineData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( fineLevel );

      const auto vertexCoarseData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( coarseLevel );
      const auto edgeCoarseData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( coarseLevel );

      if ( updateType == Replace )
      {
         // we need to set the face ghost-layers to zero explicitly since this is not necessarily done by interpolation
         for ( const auto& it : vertexdof::macroface::Iterator( fineLevel, 0 ) )
         {
            vertexFineData[vertexdof::macroface::index( fineLevel, it.x(), it.y() )] = real_c( 0 );
         }

         // For some reason the Intel compiler cannot create code if an iterator is used for the edge unknowns.
         // Therefore we use a plain for loop here.
         //
         // See issue #94.
         //
         const uint_t edgedofFieldSize = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getSize( fineLevel );
         for ( uint_t i = 0; i < edgedofFieldSize; i++ )
         {
            edgeFineData[i] = real_c( 0 );
         }
      }
      else if ( updateType == Add )
      {
         // we only set the ghost layers to zero, but not the inner unknowns
         for ( const auto& it : vertexdof::macroface::Iterator( fineLevel ) )
         {
            if ( vertexdof::macroface::isVertexOnBoundary( fineLevel, it ) )
            {
               vertexFineData[vertexdof::macroface::index( fineLevel, it.x(), it.y() )] = real_c( 0 );
            }
         }

         for ( const auto & it : edgedof::macroface::Iterator( fineLevel ) )
         {
            if ( edgedof::macroface::isHorizontalEdgeOnBoundary( fineLevel, it ) )
            {
               edgeFineData[edgedof::macroface::index( fineLevel, it.x(), it.y(), edgedof::EdgeDoFOrientation::X )] = real_c( 0 );
            }

            if ( edgedof::macroface::isDiagonalEdgeOnBoundary( fineLevel, it ) )
            {
               edgeFineData[edgedof::macroface::index( fineLevel, it.x(), it.y(), edgedof::EdgeDoFOrientation::XY )] = real_c( 0 );
            }

            if ( edgedof::macroface::isVerticalEdgeOnBoundary( fineLevel, it ) )
            {
               edgeFineData[edgedof::macroface::index( fineLevel, it.x(), it.y(), edgedof::EdgeDoFOrientation::Y )] = real_c( 0 );
            }
         }
      }
      else
      {
         WALBERLA_ABORT( "Invalid update type in prolongation." );
      }


      const double numNeighborFacesEdge0 =
          static_cast< double >( storage->getEdge( face->neighborEdges().at( 0 ) )->getNumNeighborFaces() );
      const double numNeighborFacesEdge1 =
          static_cast< double >( storage->getEdge( face->neighborEdges().at( 1 ) )->getNumNeighborFaces() );
      const double numNeighborFacesEdge2 =
          static_cast< double >( storage->getEdge( face->neighborEdges().at( 2 ) )->getNumNeighborFaces() );
      const double numNeighborFacesVertex0 =
          static_cast< double >( storage->getVertex( face->neighborVertices().at( 0 ) )->getNumNeighborFaces() );
      const double numNeighborFacesVertex1 =
          static_cast< double >( storage->getVertex( face->neighborVertices().at( 1 ) )->getNumNeighborFaces() );
      const double numNeighborFacesVertex2 =
          static_cast< double >( storage->getVertex( face->neighborVertices().at( 2 ) )->getNumNeighborFaces() );

      typedef edgedof::EdgeDoFOrientation eo;
      std::map< eo, uint_t >              firstIdxFine;
      std::map< eo, uint_t >              firstIdxCoarse;
      for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
      {
         firstIdxFine[e]   = edgedof::macroface::index( fineLevel, 0, 0, e );
         firstIdxCoarse[e] = edgedof::macroface::index( coarseLevel, 0, 0, e );
      }

      P2::macroface::generated::prolongate_2D_macroface_P2_push_from_vertexdofs( &edgeFineData[firstIdxFine[eo::X]],
                                                                                 &edgeFineData[firstIdxFine[eo::XY]],
                                                                                 &edgeFineData[firstIdxFine[eo::Y]],
                                                                                 vertexCoarseData,
                                                                                 vertexFineData,
                                                                                 static_cast< int32_t >( coarseLevel ),
                                                                                 numNeighborFacesEdge0,
                                                                                 numNeighborFacesEdge1,
                                                                                 numNeighborFacesEdge2,
                                                                                 numNeighborFacesVertex0,
                                                                                 numNeighborFacesVertex1,
                                                                                 numNeighborFacesVertex2 );

      if ( coarseLevel == 0 )
      {
         P2::macroface::generated::prolongate_2D_macroface_P2_push_from_edgedofs_level_0_to_1(
             &edgeCoarseData[firstIdxCoarse[eo::X]],
             &edgeCoarseData[firstIdxCoarse[eo::XY]],
             &edgeCoarseData[firstIdxCoarse[eo::Y]],
             &edgeFineData[firstIdxFine[eo::X]],
             &edgeFineData[firstIdxFine[eo::XY]],
             &edgeFineData[firstIdxFine[eo::Y]],
             vertexFineData,
             static_cast< int32_t >( coarseLevel ),
             numNeighborFacesEdge0,
             numNeighborFacesEdge1,
             numNeighborFacesEdge2 );
      }
      else
      {
         P2::macroface::generated::prolongate_2D_macroface_P2_push_from_edgedofs( &edgeCoarseData[firstIdxCoarse[eo::X]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::XY]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::Y]],
                                                                                  &edgeFineData[firstIdxFine[eo::X]],
                                                                                  &edgeFineData[firstIdxFine[eo::XY]],
                                                                                  &edgeFineData[firstIdxFine[eo::Y]],
                                                                                  vertexFineData,
                                                                                  static_cast< int32_t >( coarseLevel ),
                                                                                  numNeighborFacesEdge0,
                                                                                  numNeighborFacesEdge1,
                                                                                  numNeighborFacesEdge2 );
      }
   }

   function.getVertexDoFFunction().communicateAdditively< Face, Edge >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );
   function.getVertexDoFFunction().communicateAdditively< Face, Vertex >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );

   function.getEdgeDoFFunction().communicateAdditively< Face, Edge >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );
}

void P2toP2QuadraticProlongation::prolongateAdditively3D( const P2Function< real_t >& function,
                                                          const uint_t&               sourceLevel,
                                                          const DoFType&              flag,
                                                          const UpdateType &          updateType ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
   const DoFType excludeFlag = ( flag ^ All );

   const auto storage = function.getStorage();

   const uint_t fineLevel   = sourceLevel + 1;
   const uint_t coarseLevel = sourceLevel;

   function.communicate< Vertex, Edge >( coarseLevel );
   function.communicate< Edge, Face >( coarseLevel );
   function.communicate< Face, Cell >( coarseLevel );

   for ( const auto& cellIt : function.getStorage()->getCells() )
   {
      const auto cell = cellIt.second;

      auto vertexFineData = cell->getData( function.getVertexDoFFunction().getCellDataID() )->getPointer( fineLevel );
      auto edgeFineData   = cell->getData( function.getEdgeDoFFunction().getCellDataID() )->getPointer( fineLevel );

      const auto vertexCoarseData = cell->getData( function.getVertexDoFFunction().getCellDataID() )->getPointer( coarseLevel );
      const auto edgeCoarseData   = cell->getData( function.getEdgeDoFFunction().getCellDataID() )->getPointer( coarseLevel );

      if ( updateType == Replace )
      {
         // we need to set the face ghost-layers to zero explicitly since this is not necessarily done by interpolation
         for ( const auto& it : vertexdof::macrocell::Iterator( fineLevel, 0 ) )
         {
            vertexFineData[vertexdof::macrocell::index( fineLevel, it.x(), it.y(), it.z() )] = real_c( 0 );
         }

         // For some reason the Intel compiler cannot create code if an iterator is used for the edge unknowns.
         // Therefore we use a plain for loop here.
         //
         // See issue #94.
         //
         const uint_t edgedofFieldSize = cell->getData( function.getEdgeDoFFunction().getCellDataID() )->getSize( fineLevel );
         for ( uint_t i = 0; i < edgedofFieldSize; i++ )
         {
            edgeFineData[i] = real_c( 0 );
         }
      }
      else if ( updateType == Add )
      {
         // we only set the ghost layers to zero, but not the inner unknowns
         for ( const auto& it : vertexdof::macrocell::Iterator( fineLevel, 0 ) )
         {
            if ( !vertexdof::macrocell::isOnCellFace( it, fineLevel ).empty() )
            {
               vertexFineData[vertexdof::macrocell::index( fineLevel, it.x(), it.y(), it.z() )] = real_c( 0 );
            }
         }

         for ( auto orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
         {
            for ( const auto & it : edgedof::macrocell::Iterator( fineLevel ) )
            {
               if ( !edgedof::macrocell::isInnerEdgeDoF( fineLevel, it, orientation ) )
               {
                  edgeFineData[edgedof::macrocell::index( fineLevel, it.x(), it.y(), it.z(), orientation )] = real_c( 0 );
               }
            }
         }

         // no xyz edges lie on the boundary by definition
      }
      else
      {
         WALBERLA_ABORT( "Invalid update type in prolongation." );
      }


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

      typedef edgedof::EdgeDoFOrientation eo;
      std::map< eo, uint_t >              firstIdxFine;
      std::map< eo, uint_t >              firstIdxCoarse;
      for ( auto e : edgedof::allEdgeDoFOrientations )
      {
         firstIdxFine[e]   = edgedof::macrocell::index( fineLevel, 0, 0, 0, e );
         firstIdxCoarse[e] = edgedof::macrocell::index( coarseLevel, 0, 0, 0, e );
      }

      P2::macrocell::generated::prolongate_3D_macrocell_P2_push_from_vertexdofs( &edgeFineData[firstIdxFine[eo::X]],
                                                                                 &edgeFineData[firstIdxFine[eo::XY]],
                                                                                 &edgeFineData[firstIdxFine[eo::XYZ]],
                                                                                 &edgeFineData[firstIdxFine[eo::XZ]],
                                                                                 &edgeFineData[firstIdxFine[eo::Y]],
                                                                                 &edgeFineData[firstIdxFine[eo::YZ]],
                                                                                 &edgeFineData[firstIdxFine[eo::Z]],
                                                                                 vertexCoarseData,
                                                                                 vertexFineData,
                                                                                 static_cast< int32_t >( coarseLevel ),
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

      if ( coarseLevel == 0 )
      {
         P2::macrocell::generated::prolongate_3D_macrocell_P2_push_from_edgedofs_level_0_to_1(
             &edgeCoarseData[firstIdxCoarse[eo::X]],
             &edgeCoarseData[firstIdxCoarse[eo::XY]],
             &edgeCoarseData[firstIdxCoarse[eo::XZ]],
             &edgeCoarseData[firstIdxCoarse[eo::Y]],
             &edgeCoarseData[firstIdxCoarse[eo::YZ]],
             &edgeCoarseData[firstIdxCoarse[eo::Z]],
             &edgeFineData[firstIdxFine[eo::X]],
             &edgeFineData[firstIdxFine[eo::XY]],
             &edgeFineData[firstIdxFine[eo::XYZ]],
             &edgeFineData[firstIdxFine[eo::XZ]],
             &edgeFineData[firstIdxFine[eo::Y]],
             &edgeFineData[firstIdxFine[eo::YZ]],
             &edgeFineData[firstIdxFine[eo::Z]],
             vertexFineData,
             static_cast< int32_t >( coarseLevel ),
             numNeighborCellsEdge0,
             numNeighborCellsEdge1,
             numNeighborCellsEdge2,
             numNeighborCellsEdge3,
             numNeighborCellsEdge4,
             numNeighborCellsEdge5,
             numNeighborCellsFace0,
             numNeighborCellsFace1,
             numNeighborCellsFace2,
             numNeighborCellsFace3 );
      }
      else
      {
         P2::macrocell::generated::prolongate_3D_macrocell_P2_push_from_edgedofs( &edgeCoarseData[firstIdxCoarse[eo::X]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::XY]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::XYZ]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::XZ]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::Y]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::YZ]],
                                                                                  &edgeCoarseData[firstIdxCoarse[eo::Z]],
                                                                                  &edgeFineData[firstIdxFine[eo::X]],
                                                                                  &edgeFineData[firstIdxFine[eo::XY]],
                                                                                  &edgeFineData[firstIdxFine[eo::XYZ]],
                                                                                  &edgeFineData[firstIdxFine[eo::XZ]],
                                                                                  &edgeFineData[firstIdxFine[eo::Y]],
                                                                                  &edgeFineData[firstIdxFine[eo::YZ]],
                                                                                  &edgeFineData[firstIdxFine[eo::Z]],
                                                                                  vertexFineData,
                                                                                  static_cast< int32_t >( coarseLevel ),
                                                                                  numNeighborCellsEdge0,
                                                                                  numNeighborCellsEdge1,
                                                                                  numNeighborCellsEdge2,
                                                                                  numNeighborCellsEdge3,
                                                                                  numNeighborCellsEdge4,
                                                                                  numNeighborCellsEdge5,
                                                                                  numNeighborCellsFace0,
                                                                                  numNeighborCellsFace1,
                                                                                  numNeighborCellsFace2,
                                                                                  numNeighborCellsFace3 );
      }
   }

   function.getVertexDoFFunction().communicateAdditively< Cell, Face >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );
   function.getVertexDoFFunction().communicateAdditively< Cell, Edge >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );
   function.getVertexDoFFunction().communicateAdditively< Cell, Vertex >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );

   function.getEdgeDoFFunction().communicateAdditively< Cell, Face >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );
   function.getEdgeDoFFunction().communicateAdditively< Cell, Edge >( fineLevel, excludeFlag, *function.getStorage(), updateType == Replace );
}

void P2toP2QuadraticProlongation::prolongateStandard( const P2Function< real_t >& function,
                                                      const uint_t&               sourceLevel,
                                                      const DoFType&              flag ) const
{
   const auto& vertexDoFFunction = function.getVertexDoFFunction();
   const auto& edgeDoFFunction   = function.getEdgeDoFFunction();

   edgeDoFFunction.template communicate< Vertex, Edge >( sourceLevel );
   edgeDoFFunction.template communicate< Edge, Face >( sourceLevel );

   vertexDoFFunction.template communicate< Vertex, Edge >( sourceLevel );
   vertexDoFFunction.template communicate< Edge, Face >( sourceLevel );

   for ( const auto& it : function.getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::prolongate< real_t >(
             sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
      }
   }

   for ( const auto& it : function.getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::prolongate< real_t >(
             sourceLevel, edge, vertexDoFFunction.getEdgeDataID(), edgeDoFFunction.getEdgeDataID() );
      }
   }

   for ( const auto& it : function.getStorage()->getVertices() )
   {
      const Vertex& vertex = *it.second;

      const DoFType vertexBC = function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::prolongate< real_t >(
             sourceLevel, vertex, vertexDoFFunction.getVertexDataID(), edgeDoFFunction.getVertexDataID() );
      }
   }
}
} // namespace hyteg
