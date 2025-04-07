/*
 * Copyright (c) 2024 Marcus Mohr.
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

/// \file Utility for exporting a refined mesh from HyTeG in MSH4.1 format
#pragma once

#include <fstream>
#include <string>

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"

namespace hyteg::gmsh {

/// Helper function for exportRefinedMesh
void writeMeshSection( std::ofstream& meshFile )
{
   meshFile << "$MeshFormat\n4.1 0 8\n$EndMeshFormat" << std::endl;
}

/// Helper function for exportRefinedMesh
void writeNodesSection( std::ofstream&                       meshFile,
                        std::shared_ptr< PrimitiveStorage >& storage,
                        uint_t                               exportLevel,
                        P1Function< idx_t >&                 numerator )
{
   // section start marker
   meshFile << "$Nodes" << std::endl;

   // number of entities corresponds to number of primitives in storage
   uint_t numEntities = storage->getNumberOfGlobalPrimitives();

   // number of nodes corresponds to number of DoFs for a P1 function on the level,
   // we consecutively number them from 1 to how many there are
   uint_t numNodes   = numberOfGlobalDoFs< VertexDoFFunctionTag >( *storage, exportLevel );
   uint_t minNodeTag = 1u;
   uint_t maxNodeTag = numNodes;

   meshFile << numEntities << " " << numNodes << " " << minNodeTag << " " << maxNodeTag << std::endl;

   uint_t nodeTag = minNodeTag;

   // export nodes on macro-vertices
   uint_t vertexTag = 1u;

   for ( const auto& vertexIter : storage->getVertices() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing vertex with tag: " << vertexTag );

      auto vertexPtr = vertexIter.second;

      // obtain micro-vertex coordinates in physical domain
      const Point3D coordinate = vertexPtr->getCoordinates();
      Point3D       blendedCoords;
      vertexPtr->getGeometryMap()->evalF( coordinate, blendedCoords );

      // obtain global node index
      auto numeratorData = vertexPtr->getData( numerator.getVertexDataID() )->getPointer( exportLevel );
      nodeTag            = numeratorData[0];

      // dimension -- nodeTag -- parametric -- number of nodes in block
      meshFile << "0 " << vertexTag << " 0 1\n";
      meshFile << nodeTag << "\n";
      meshFile << blendedCoords[0] << " " << blendedCoords[1] << " " << blendedCoords[2] << std::endl;
      ++vertexTag;
   }

   // export nodes on macro-edges
   uint_t edgeTag = 1u;
   for ( const auto& edgeIter : storage->getEdges() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing edge with tag: " << edgeTag );

      auto edgePtr = edgeIter.second;

      uint_t nodesPerEdge = levelinfo::num_microvertices_per_edge( exportLevel ) - 2;

      // dimension -- edgeTag -- parametric -- number of nodes in block
      meshFile << "1 " << edgeTag << " 0 " << nodesPerEdge << "\n";

      if ( nodesPerEdge > 0 )
      {
         for ( const auto& it : vertexdof::macroedge::Iterator( exportLevel, 1 ) )
         {
            // obtain global node index
            const uint_t idx           = vertexdof::macroedge::indexFromVertex( exportLevel, it.x(), stencilDirection::VERTEX_C );
            auto         numeratorData = edgePtr->getData( numerator.getEdgeDataID() )->getPointer( exportLevel );
            nodeTag                    = numeratorData[idx];

            // write out data
            meshFile << nodeTag << "\n";
         }

         for ( const auto& it : vertexdof::macroedge::Iterator( exportLevel, 1 ) )
         {
            // obtain micro-vertex coordinates in physical domain
            const Point3D coordinate = vertexdof::macroedge::coordinateFromIndex( exportLevel, *edgePtr, it );
            Point3D       blendedCoords;
            edgePtr->getGeometryMap()->evalF( coordinate, blendedCoords );

            // write out data
            meshFile << blendedCoords[0] << " " << blendedCoords[1] << " " << blendedCoords[2] << std::endl;
         }
      }

      ++edgeTag;
   }

   // export nodes on macro-faces
   uint_t faceTag = 1u;
   for ( const auto& faceIter : storage->getFaces() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing face with tag: " << faceTag );

      auto facePtr = faceIter.second;

      uint_t nodesPerFace = numberOfInnerDoFs< P1FunctionTag, Face >( exportLevel );

      // dimension -- faceTag -- parametric -- number of nodes in block
      meshFile << "2 " << faceTag << " 0 " << nodesPerFace << "\n";

      if ( nodesPerFace > 0 )
      {
         for ( const auto& it : vertexdof::macroface::Iterator( exportLevel, 1 ) )
         {
            // obtain global node index
            const uint_t idx = vertexdof::macroface::indexFromVertex( exportLevel, it.x(), it.y(), stencilDirection::VERTEX_C );
            auto         numeratorData = facePtr->getData( numerator.getFaceDataID() )->getPointer( exportLevel );
            nodeTag                    = numeratorData[idx];

            // write out data
            meshFile << nodeTag << "\n";
         }

         for ( const auto& it : vertexdof::macroface::Iterator( exportLevel, 1 ) )
         {
            // obtain micro-vertex coordinates in physical domain
            const Point3D coordinate = vertexdof::macroface::coordinateFromIndex( exportLevel, *facePtr, it );
            Point3D       blendedCoords;
            facePtr->getGeometryMap()->evalF( coordinate, blendedCoords );

            // write out data
            meshFile << blendedCoords[0] << " " << blendedCoords[1] << " " << blendedCoords[2] << std::endl;
         }
      }

      ++faceTag;
   }

   // export nodes on macro-cells
   uint_t cellTag = 1u;
   for ( const auto& cellIter : storage->getCells() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing cell with tag: " << cellTag );

      auto cellPtr = cellIter.second;

      uint_t nodesPerCell = numberOfInnerDoFs< P1FunctionTag, Cell >( exportLevel );

      // dimension -- cellTag -- parametric -- number of nodes in block
      meshFile << "3 " << cellTag << " 0 " << nodesPerCell << "\n";

      if ( nodesPerCell > 0 )
      {
         for ( const auto& it : vertexdof::macrocell::Iterator( exportLevel, 1 ) )
         {
            // obtain global node index
            const uint_t idx =
                vertexdof::macrocell::indexFromVertex( exportLevel, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
            auto numeratorData = cellPtr->getData( numerator.getCellDataID() )->getPointer( exportLevel );
            nodeTag            = numeratorData[idx];

            // write out data
            meshFile << nodeTag << "\n";
         }

         for ( const auto& it : vertexdof::macrocell::Iterator( exportLevel, 1 ) )
         {
            // obtain micro-vertex coordinates in physical domain
            const Point3D coordinate = vertexdof::macrocell::coordinateFromIndex( exportLevel, *cellPtr, it );
            Point3D       blendedCoords;
            cellPtr->getGeometryMap()->evalF( coordinate, blendedCoords );

            // write out data
            meshFile << blendedCoords[0] << " " << blendedCoords[1] << " " << blendedCoords[2] << std::endl;
         }
      }
      ++cellTag;
   }

   // section stop marker
   meshFile << "$EndNodes" << std::endl;
}

/// Helper function for exportRefinedMesh
void writeElementsSection( std::ofstream&                       meshFile,
                           std::shared_ptr< PrimitiveStorage >& storage,
                           uint_t                               exportLevel,
                           P1Function< idx_t >&                 numerator )
{
   // section start marker
   meshFile << "$Elements" << std::endl;

   // determine number of elements for Gmsh, i.e. number of micro-elements of the same dimension
   // as the associated primitive
   uint_t linesPerEdge     = levelinfo::num_microedges_per_edge( exportLevel );
   uint_t trianglesPerFace = levelinfo::num_microfaces_per_face( exportLevel );
   uint_t tetsPerCell      = levelinfo::num_microcells_per_cell( exportLevel );

   uint_t numElements{ 0 };

   numElements += storage->getNumberOfLocalVertices();
   numElements += storage->getNumberOfLocalEdges() * linesPerEdge;
   numElements += storage->getNumberOfLocalFaces() * trianglesPerFace;
   numElements += storage->getNumberOfLocalCells() * tetsPerCell;

   WALBERLA_LOG_DETAIL_ON_ROOT( "Need to write " << numElements << " elements into MSH file." );

   // numEntityBlocks -- numElements -- minElementTag -- maxElementTag
   meshFile << storage->getNumberOfLocalPrimitives() << " " << numElements << " 1 " << numElements << std::endl;

   // write node elements from macro-vertices
   uint_t vertexTag  = 1u;
   uint_t elementTag = 1u;

   for ( const auto& vertexIter : storage->getVertices() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing vertex with tag: " << vertexTag );

      auto vertexPtr = vertexIter.second;

      // obtain global node index
      auto  numeratorData = vertexPtr->getData( numerator.getVertexDataID() )->getPointer( exportLevel );
      idx_t nodeTag       = numeratorData[0];

      // entitiyDim -- entityTag -- elementType -- numElementsInBlock
      meshFile << "0 " << vertexTag << " 15 1\n";

      // elementTag -- nodeTag
      meshFile << elementTag << " " << nodeTag << "\n";
      ++vertexTag;
      ++elementTag;
   }

   // write line elements from macro-edges
   uint_t edgeTag = 1u;

   for ( const auto& edgeIter : storage->getEdges() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing edge with tag: " << edgeTag );

      auto edgePtr = edgeIter.second;

      // entitiyDim -- entityTag -- elementType -- numElementsInBlock
      meshFile << "1 " << edgeTag << " 1 " << linesPerEdge << "\n";

      if ( linesPerEdge > 0 )
      {
         for ( const auto& it : edgedof::macroedge::Iterator( exportLevel ) )
         {
            const uint_t leftIdx  = vertexdof::macroedge::indexFromVertex( exportLevel, it.x(), stencilDirection::VERTEX_C );
            const uint_t rightIdx = vertexdof::macroedge::indexFromVertex( exportLevel, it.x(), stencilDirection::VERTEX_E );

            auto numeratorData = edgePtr->getData( numerator.getEdgeDataID() )->getPointer( exportLevel );

            uint_t leftTag  = numeratorData[leftIdx];
            uint_t rightTag = numeratorData[rightIdx];

            // elementTag -- nodeTags
            meshFile << elementTag << " " << leftTag << " " << rightTag << "\n";

            ++elementTag;
         }
      }

      ++edgeTag;
   }

   // write triangle elements from macro-faces
   uint_t faceTag = 1u;

   for ( const auto& faceIter : storage->getFaces() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing face with tag: " << faceTag );

      auto facePtr = faceIter.second;

      // entitiyDim -- entityTag -- elementType -- numElementsInBlock
      meshFile << "2 " << faceTag << " 2 " << trianglesPerFace << "\n";

      if ( trianglesPerFace > 0 )
      {
         // loop over micro-faces
         for ( const auto& faceType : facedof::allFaceTypes )
         {
            for ( const auto& microFace : facedof::macroface::Iterator( exportLevel, faceType, 0 ) )
            {
               std::array< uint_t, 3 > vertexDoFIndices;
               vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, exportLevel, vertexDoFIndices );

               auto numeratorData = facePtr->getData( numerator.getFaceDataID() )->getPointer( exportLevel );

               idx_t node0 = numeratorData[vertexDoFIndices[0]];
               idx_t node1 = numeratorData[vertexDoFIndices[1]];
               idx_t node2 = numeratorData[vertexDoFIndices[2]];

               // elementTag -- nodeTags
               meshFile << elementTag << " " << node0 << " " << node1 << " " << node2 << "\n";
               ++elementTag;
            }
         }
      }
      ++faceTag;
   }

   // write tetrahedron elements from macro-cells
   uint_t cellTag = 1u;

   for ( const auto& cellIter : storage->getCells() )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( "-> Processing cell with tag: " << cellTag );

      auto cellPtr = cellIter.second;

      // entitiyDim -- entityTag -- elementType -- numElementsInBlock
      meshFile << "3 " << cellTag << " 4 " << tetsPerCell << "\n";

      if ( tetsPerCell > 0 )
      {
         // loop over micro-cells
         for ( const auto& cellType : celldof::allCellTypes )
         {
            for ( const auto& microCell : celldof::macrocell::Iterator( exportLevel, cellType, 0 ) )
            {
               std::array< uint_t, 4 > vertexDoFIndices;
               vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cellType, exportLevel, vertexDoFIndices );

               auto numeratorData = cellPtr->getData( numerator.getCellDataID() )->getPointer( exportLevel );

               idx_t node0 = numeratorData[vertexDoFIndices[0]];
               idx_t node1 = numeratorData[vertexDoFIndices[1]];
               idx_t node2 = numeratorData[vertexDoFIndices[2]];
               idx_t node3 = numeratorData[vertexDoFIndices[3]];

               // elementTag -- nodeTags
               meshFile << elementTag << " " << node0 << " " << node1 << " " << node2 << " " << node3 << "\n";
               ++elementTag;
            }
         }
      }

      ++cellTag;
   }

   // section stop marker
   meshFile << "$EndElements" << std::endl;
}

/// Free-function for exporting a refined HyTeG mesh to a file in MSH4.1 format
///
/// \note This function only works in a sequential setting. Currently it cannot be used to export a mesh
///       distributed over multiple MPI processes.
///
/// \param storage         A shared pointer to a PrimitiveStorage object containing the base mesh
/// \param exportLevel     Refinement level for which we want to export the mesh
/// \param exportFileName  Name of the MSH4.1 file to generate
void exportRefinedMesh( std::shared_ptr< PrimitiveStorage >& storage, uint_t exportLevel, std::string exportFileName )
{
   std::ofstream meshFile( exportFileName );

   if ( !meshFile )
   {
      WALBERLA_ABORT( "gmsh::exportRefinedMesh(): Failed to open '" << exportFileName << "'" );
   }

#ifdef HYTEG_BUILD_WITH_MPI
   if ( walberla::mpi::MPIManager::instance()->numProcesses() > 1 )
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
          "gmsh::exportRefinedMesh(): is not designed to work with multiple MPI processes! Proceed at your own risk!" );
   }
#endif

   // The DoFs of a P1Function correspond to the micro-vertices of the mesh.
   // Hence, we can use an enumerator to index the micro-vertices.
   P1Function< idx_t > numerator( "Gmsh numerator", storage, exportLevel, exportLevel );
   idx_t               offset = static_cast< idx_t >( 1 );
   numerator.enumerate( exportLevel, offset );

   WALBERLA_LOG_INFO_ON_ROOT( "Exporting mesh of refinement level " << exportLevel << " to file '" << exportFileName << "'" );

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Going to write $Mesh section." );
   writeMeshSection( meshFile );

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Going to write $Nodes section." );
   writeNodesSection( meshFile, storage, exportLevel, numerator );

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Going to write $Elements section." );
   writeElementsSection( meshFile, storage, exportLevel, numerator );

   meshFile.close();
}

} // namespace hyteg::gmsh
