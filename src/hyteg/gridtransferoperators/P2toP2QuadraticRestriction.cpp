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
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"

#include "hyteg/gridtransferoperators/generatedKernels/restrict_2D_macroface_P2_update_edgedofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/restrict_2D_macroface_P2_update_edgedofs_level_0_to_1.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/restrict_2D_macroface_P2_update_vertexdofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/restrict_3D_macrocell_P2_update_vertexdofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/restrict_3D_macrocell_P2_update_edgedofs.hpp"
#include "hyteg/gridtransferoperators/generatedKernels/restrict_3D_macrocell_P2_update_edgedofs_level_1_to_0.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p2functionspace/P2Multigrid.hpp"

namespace hyteg {

using indexing::Index;

void P2toP2QuadraticRestriction::restrictAdditively( const P2Function< real_t >& function,
                                                     const uint_t&               sourceLevel,
                                                     const DoFType&              flag ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
   const DoFType excludeFlag = ( flag ^ All );

   const auto storage = function.getStorage();

   const uint_t fineLevel   = sourceLevel;
   const uint_t coarseLevel = sourceLevel - 1;

   function.communicate< Vertex, Edge >( fineLevel );
   function.communicate< Edge, Face >( fineLevel );

   for ( const auto& faceIt : function.getStorage()->getFaces() )
   {
      const auto face = faceIt.second;

      const auto vertexFineData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( fineLevel );
      const auto edgeFineData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( fineLevel );

      auto vertexCoarseData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( coarseLevel );
      auto edgeCoarseData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( coarseLevel );

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

      typedef edgedof::EdgeDoFOrientation eo;
      std::map< eo, uint_t >              firstIdxFine;
      std::map< eo, uint_t >              firstIdxCoarse;
      for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
      {
         firstIdxFine[e]   = edgedof::macroface::index( fineLevel, 0, 0, e );
         firstIdxCoarse[e] = edgedof::macroface::index( coarseLevel, 0, 0, e );
      }

      P2::macroface::generated::restrict_2D_macroface_P2_update_vertexdofs( &edgeFineData[firstIdxFine[eo::X]],
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
         P2::macroface::generated::restrict_2D_macroface_P2_update_edgedofs_level_0_to_1( &edgeCoarseData[firstIdxCoarse[eo::X]],
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
         P2::macroface::generated::restrict_2D_macroface_P2_update_edgedofs( &edgeCoarseData[firstIdxCoarse[eo::X]],
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

   function.getVertexDoFFunction().communicateAdditively< Face, Edge >( coarseLevel, excludeFlag, *function.getStorage() );
   function.getVertexDoFFunction().communicateAdditively< Face, Vertex >( coarseLevel, excludeFlag, *function.getStorage() );

   function.getEdgeDoFFunction().communicateAdditively< Face, Edge >( coarseLevel, excludeFlag, *function.getStorage() );
}

void P2toP2QuadraticRestriction::restrictAdditively3D( const P2Function< real_t >& function,
                                                       const uint_t&               sourceLevel,
                                                       const DoFType&              flag ) const
{
   /// XOR flag with all to get the DoFTypes that should be excluded
   const DoFType excludeFlag = ( flag ^ All );
   const auto    storage     = function.getStorage();

   const uint_t fineLevel   = sourceLevel;
   const uint_t coarseLevel = sourceLevel - 1;

   function.communicate< Vertex, Edge >( fineLevel );
   function.communicate< Edge, Face >( fineLevel );
   function.communicate< Face, Cell >( fineLevel );

   for ( const auto& cellIt : function.getStorage()->getCells() )
   {
      const auto cell = cellIt.second;

      const auto vertexFineData = cell->getData( function.getVertexDoFFunction().getCellDataID() )->getPointer( fineLevel );
      const auto edgeFineData   = cell->getData( function.getEdgeDoFFunction().getCellDataID() )->getPointer( fineLevel );

      auto vertexCoarseData = cell->getData( function.getVertexDoFFunction().getCellDataID() )->getPointer( coarseLevel );
      auto edgeCoarseData   = cell->getData( function.getEdgeDoFFunction().getCellDataID() )->getPointer( coarseLevel );

      const auto numNeighborCellsFace0 =
          static_cast< double >( storage->getFace( cell->neighborFaces().at( 0 ) )->getNumNeighborCells() );
      const auto numNeighborCellsFace1 =
          static_cast< double >( storage->getFace( cell->neighborFaces().at( 1 ) )->getNumNeighborCells() );
      const auto numNeighborCellsFace2 =
          static_cast< double >( storage->getFace( cell->neighborFaces().at( 2 ) )->getNumNeighborCells() );
      const auto numNeighborCellsFace3 =
          static_cast< double >( storage->getFace( cell->neighborFaces().at( 3 ) )->getNumNeighborCells() );

      const auto numNeighborCellsEdge0 =
          static_cast< double >( storage->getEdge( cell->neighborEdges().at( 0 ) )->getNumNeighborCells() );
      const auto numNeighborCellsEdge1 =
          static_cast< double >( storage->getEdge( cell->neighborEdges().at( 1 ) )->getNumNeighborCells() );
      const auto numNeighborCellsEdge2 =
          static_cast< double >( storage->getEdge( cell->neighborEdges().at( 2 ) )->getNumNeighborCells() );
      const auto numNeighborCellsEdge3 =
          static_cast< double >( storage->getEdge( cell->neighborEdges().at( 3 ) )->getNumNeighborCells() );
      const auto numNeighborCellsEdge4 =
          static_cast< double >( storage->getEdge( cell->neighborEdges().at( 4 ) )->getNumNeighborCells() );
      const auto numNeighborCellsEdge5 =
          static_cast< double >( storage->getEdge( cell->neighborEdges().at( 5 ) )->getNumNeighborCells() );

      const auto numNeighborCellsVertex0 =
          static_cast< double >( storage->getVertex( cell->neighborVertices().at( 0 ) )->getNumNeighborCells() );
      const auto numNeighborCellsVertex1 =
          static_cast< double >( storage->getVertex( cell->neighborVertices().at( 1 ) )->getNumNeighborCells() );
      const auto numNeighborCellsVertex2 =
          static_cast< double >( storage->getVertex( cell->neighborVertices().at( 2 ) )->getNumNeighborCells() );
      const auto numNeighborCellsVertex3 =
          static_cast< double >( storage->getVertex( cell->neighborVertices().at( 3 ) )->getNumNeighborCells() );

      typedef edgedof::EdgeDoFOrientation eo;
      std::map< eo, uint_t >              firstIdxFine;
      std::map< eo, uint_t >              firstIdxCoarse;
      for ( auto e : edgedof::allEdgeDoFOrientations )
      {
         firstIdxFine[e]   = edgedof::macrocell::index( fineLevel, 0, 0, 0, e );
         firstIdxCoarse[e] = edgedof::macrocell::index( coarseLevel, 0, 0, 0, e );
      }

      P2::macrocell::generated::restrict_3D_macrocell_P2_update_vertexdofs( &edgeFineData[firstIdxFine[eo::X]],
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
         P2::macrocell::generated::restrict_3D_macrocell_P2_update_edgedofs_level_1_to_0( &edgeCoarseData[firstIdxCoarse[eo::X]],
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
         P2::macrocell::generated::restrict_3D_macrocell_P2_update_edgedofs( &edgeCoarseData[firstIdxCoarse[eo::X]],
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

   function.getVertexDoFFunction().communicateAdditively< Cell, Face >( coarseLevel, excludeFlag, *function.getStorage() );
   function.getVertexDoFFunction().communicateAdditively< Cell, Edge >( coarseLevel, excludeFlag, *function.getStorage() );
   function.getVertexDoFFunction().communicateAdditively< Cell, Vertex >( coarseLevel, excludeFlag, *function.getStorage() );

   function.getEdgeDoFFunction().communicateAdditively< Cell, Face >( coarseLevel, excludeFlag, *function.getStorage() );
   function.getEdgeDoFFunction().communicateAdditively< Cell, Edge >( coarseLevel, excludeFlag, *function.getStorage() );
}

void P2toP2QuadraticRestriction::restrictWithPostCommunication( const hyteg::P2Function< walberla::real_t >& function,
                                                                const uint_t&                                sourceLevel,
                                                                const hyteg::DoFType&                        flag ) const
{
   const auto  storage           = function.getStorage();
   const auto& vertexDoFFunction = function.getVertexDoFFunction();
   const auto& edgeDoFFunction   = function.getEdgeDoFFunction();
   const auto  boundaryCondition = function.getBoundaryCondition();

   edgeDoFFunction.template communicate< Vertex, Edge >( sourceLevel );
   edgeDoFFunction.template communicate< Edge, Face >( sourceLevel );

   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::restrict< real_t >(
             sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
      }
   }

   /// sync the vertex dofs which contain the missing edge dofs
   edgeDoFFunction.template communicate< Face, Edge >( sourceLevel );

   /// remove the temporary updates
   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::postRestrict< real_t >(
             sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
      }
   }

   for ( const auto& it : storage->getEdges() )
   {
      const Edge& edge = *it.second;

      const DoFType edgeBC = boundaryCondition.getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::restrict< real_t >(
             sourceLevel, edge, vertexDoFFunction.getEdgeDataID(), edgeDoFFunction.getEdgeDataID() );
      }
   }

   //TODO: add real vertex restrict
   for ( const auto& it : storage->getVertices() )
   {
      const Vertex& vertex = *it.second;

      const DoFType vertexBC = boundaryCondition.getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::restrictInjection< real_t >(
             sourceLevel, vertex, vertexDoFFunction.getVertexDataID(), edgeDoFFunction.getVertexDataID() );
      }
   }
}

} // namespace hyteg
