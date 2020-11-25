/*
 * Copyright (c) 2020 Marcus Mohr.
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
#include <set>
#include <sstream>

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/edgedofspace/EdgeDoFOperator.hpp"
#include "hyteg/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/indexing/CouplingCount.hpp"
// #include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p2functionspace/P2Elements3D.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

// #define SEPARATE_CELL_COUPLINGS_FOR_EDGES_AND_FACES
// #define SEPARATE_FACE_COUPLINGS_FOR_FACES
// #define SEPARATE_EDGE_COUPLINGS_FOR_EDGES

// ==================================================================================================================================

uint_t getTotalCouplings( uint_t level )
{
   uint_t numInCellCouplings = 0;
   for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) )
   {
      std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

      if ( edgedof::macrocell::isInnerXEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
      if ( edgedof::macrocell::isInnerYEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
      if ( edgedof::macrocell::isInnerZEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
      if ( edgedof::macrocell::isInnerXYEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
      if ( edgedof::macrocell::isInnerXZEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
      if ( edgedof::macrocell::isInnerYZEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

      for ( const auto& centerOrientation : innerOrientations )
      {
         for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
         {
            numInCellCouplings +=
                P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation )
                    .size();
         }
      }
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) )
   {
      WALBERLA_UNUSED( it );
      const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;
      for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
      {
         numInCellCouplings +=
             P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation ).size();
      }
   }

   return numInCellCouplings;
}

// ==================================================================================================================================

// returns an array with total coupling number, cell-to-cell, cell-to-face and cell-to-edge coupling numbers
std::array< uint_t, 3 > getCellCouplingsByType( uint_t level )
{
   uint_t numTotalCouplings     = 0;
   uint_t numCell2CellCouplings = 0;
   uint_t numCell2FaceCouplings = 0;
   uint_t numCell2EdgeCouplings = 0;

#ifdef SEPARATE_CELL_COUPLINGS_FOR_EDGES_AND_FACES
   std::array< uint_t, 4 > numCell2FacePerFace = {0};
   std::array< uint_t, 6 > numCell2EdgePerEdge = {0};
#endif

   for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) )
   {
      std::vector< edgedof::EdgeDoFOrientation > innerOrientations;

      if ( edgedof::macrocell::isInnerXEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::X );
      if ( edgedof::macrocell::isInnerYEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Y );
      if ( edgedof::macrocell::isInnerZEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::Z );
      if ( edgedof::macrocell::isInnerXYEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XY );
      if ( edgedof::macrocell::isInnerXZEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::XZ );
      if ( edgedof::macrocell::isInnerYZEdgeDoF( level, it ) )
         innerOrientations.push_back( edgedof::EdgeDoFOrientation::YZ );

      for ( const auto& centerOrientation : innerOrientations )
      {
         for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
         {
            const auto edgeDoFNeighbors =
                P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
            for ( const auto& neighbor : edgeDoFNeighbors )
            {
               const auto srcIdx = it + neighbor;
               if ( edgedof::macrocell::isInnerEdgeDoF( level, srcIdx, leafOrientation ) )
               {
                  numCell2CellCouplings++;
               }
               else if ( edgedof::macrocell::isOnCellEdges( level, srcIdx, leafOrientation ).size() == 0 )
               {
                  numCell2FaceCouplings++;

#ifdef SEPARATE_CELL_COUPLINGS_FOR_EDGES_AND_FACES
                  // count faces separately
                  auto faceSet = edgedof::macrocell::isOnCellFaces( level, srcIdx, leafOrientation );
                  for ( const auto& face : faceSet )
                  {
                     numCell2FacePerFace[face] += 1;
                  }
#endif
               }
               else
               {
#ifdef SEPARATE_CELL_COUPLINGS_FOR_EDGES_AND_FACES
                  // count edges separately
                  auto edgeSet = edgedof::macrocell::isOnCellEdges( level, srcIdx, leafOrientation );
                  for ( const auto& edge : edgeSet )
                  {
                     numCell2EdgePerEdge[edge] += 1;
                  }
#endif
                  numCell2EdgeCouplings++;
               }
               numTotalCouplings++;
            }
         }
      }
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) )
   {
      WALBERLA_UNUSED( it );
      const auto centerOrientation = edgedof::EdgeDoFOrientation::XYZ;
      for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
      {
         const auto edgeDoFNeighbors =
             P2Elements::P2Elements3D::getAllEdgeDoFNeighborsFromEdgeDoFInMacroCell( centerOrientation, leafOrientation );
         for ( const auto& neighbor : edgeDoFNeighbors )
         {
            const auto srcIdx = it + neighbor;
            if ( edgedof::macrocell::isInnerEdgeDoF( level, srcIdx, leafOrientation ) )
            {
               numCell2CellCouplings++;
            }
            else if ( edgedof::macrocell::isOnCellEdges( level, srcIdx, leafOrientation ).size() == 0 )
            {
#ifdef SEPARATE_CELL_COUPLINGS_FOR_EDGES_AND_FACES
               // count faces separately
               auto faceSet = edgedof::macrocell::isOnCellFaces( level, srcIdx, leafOrientation );
               for ( const auto& face : faceSet )
               {
                  numCell2FacePerFace[face] += 1;
               }
#endif
               numCell2FaceCouplings++;
            }
            else
            {
#ifdef SEPARATE_CELL_COUPLINGS_FOR_EDGES_AND_FACES
               // count edges separately
               auto edgeSet = edgedof::macrocell::isOnCellEdges( level, srcIdx, leafOrientation );
               for ( const auto& edge : edgeSet )
               {
                  numCell2EdgePerEdge[edge] += 1;
               }
#endif
               numCell2EdgeCouplings++;
            }
            numTotalCouplings++;
         }
      }
   }

#ifdef SEPARATE_CELL_COUPLINGS_FOR_EDGES_AND_FACES
   for ( uint_t k = 0; k < 4; k++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "level = " << level << ": numCell2FacePerFace[" << k << "] = " << numCell2FacePerFace[k] );
   }
   WALBERLA_LOG_INFO_ON_ROOT( ".................................................." );
   for ( uint_t k = 0; k < 6; k++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "level = " << level << ": numCell2EdgePerEdge[" << k << "] = " << numCell2EdgePerEdge[k] );
   }
   WALBERLA_LOG_INFO_ON_ROOT( ".................................................." );
#endif

   // return counts;
   WALBERLA_ASSERT_EQUAL( numTotalCouplings, numCell2CellCouplings + numCell2FaceCouplings + numCell2EdgeCouplings );
   return {numCell2CellCouplings, numCell2FaceCouplings, numCell2EdgeCouplings};
}

// ==================================================================================================================================

uint_t getFaceCouplingsForFace( const uint_t&                                                                       level,
                                Face&                                                                               face,
                                const PrimitiveStorage&                                                             storage,
                                const PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >& operatorId )
{
   uint_t numTotalCouplings = 0;
   // uint_t numFace2CellCouplings = 0;
   // uint_t numFace2FaceCouplings = 0;
   // uint_t numFace2EdgeCouplings = 0;

   auto opr_data = face.getData( operatorId )->getData( level );

   for ( const auto& centerIndexInFace : hyteg::edgedof::macroface::Iterator( level, 0 ) )
   {
      for ( const auto& faceCenterOrientation : edgedof::faceLocalEdgeDoFOrientations )
      {
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::X &&
              edgedof::macroface::isHorizontalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::Y &&
              edgedof::macroface::isVerticalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;
         if ( faceCenterOrientation == edgedof::EdgeDoFOrientation::XY &&
              edgedof::macroface::isDiagonalEdgeOnBoundary( level, centerIndexInFace ) )
            continue;

         for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
         {
            const Cell&  neighborCell = *( storage.getCell( face.neighborCells().at( neighborCellID ) ) );
            const uint_t localFaceID  = neighborCell.getLocalFaceID( face.getID() );

            const auto centerIndexInCell =
                edgedof::macroface::getIndexInNeighboringMacroCell( centerIndexInFace, face, neighborCellID, storage, level );
            const auto cellCenterOrientation =
                edgedof::macroface::getOrientattionInNeighboringMacroCell( faceCenterOrientation, face, neighborCellID, storage );

            for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
               for ( const auto& stencilIt : opr_data[neighborCellID][cellCenterOrientation][leafOrientation] )
               {
                  const auto stencilOffset = stencilIt.first;

                  const auto leafOrientationInFace = edgedof::macrocell::getOrientattionInNeighboringMacroFace(
                      leafOrientation, neighborCell, localFaceID, storage );

                  const auto leafIndexInCell = centerIndexInCell + stencilOffset;
                  const auto leafIndexInFace = leafOrientation == edgedof::EdgeDoFOrientation::XYZ ?
                                                   edgedof::macrocell::getIndexInNeighboringMacroFaceXYZ(
                                                       leafIndexInCell, neighborCell, localFaceID, storage, level ) :
                                                   edgedof::macrocell::getIndexInNeighboringMacroFace(
                                                       leafIndexInCell, neighborCell, localFaceID, storage, level );

                  WALBERLA_ASSERT_LESS_EQUAL( leafIndexInFace.z(), 1 )

                  uint_t leafArrayIndexInFace;
                  if ( algorithms::contains( edgedof::faceLocalEdgeDoFOrientations, leafOrientationInFace ) &&
                       leafIndexInFace.z() == 0 )
                  {
                     leafArrayIndexInFace =
                         edgedof::macroface::index( level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace );
                  }
                  else
                  {
                     leafArrayIndexInFace = edgedof::macroface::index(
                         level, leafIndexInFace.x(), leafIndexInFace.y(), leafOrientationInFace, neighborCellID );
                  }

                  numTotalCouplings++;
               }
            }
         }
      }
   }

   return numTotalCouplings++;
}

uint_t getFaceCouplingsTotal( const uint_t&                                                                       level,
                              PrimitiveStorage&                                                                   storage,
                              const PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >& operatorID )
{
   std::array< uint_t, 4 > perFaceCouplings = {0};

   uint_t k = 0;
   for ( PrimitiveID& faceID : storage.getFaceIDs() )
   {
      Face& face          = *( storage.getFace( faceID ) );
      perFaceCouplings[k] = getFaceCouplingsForFace( level, face, storage, operatorID );
      k++;
#ifdef SEPARATE_FACE_COUPLINGS_FOR_FACES
      WALBERLA_LOG_INFO_ON_ROOT( "couplings for current face " << k << " = " << perFaceCouplings[k-1] );
#else
      break;
#endif
   }
   return perFaceCouplings[0];
}

// ==================================================================================================================================

uint_t getEdgeCouplingsForEdge( const uint_t&                                                                       level,
                                const Edge&                                                                         edge,
                                const PrimitiveStorage&                                                             storage,
                                const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >& operatorId )
{

   auto opr_data = edge.getData( operatorId )->getData( level );
   uint_t count = 0;

   for ( const auto& centerIndexOnEdge : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      const edgedof::EdgeDoFOrientation edgeCenterOrientation = edgedof::EdgeDoFOrientation::X;

      for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
      {
         const Cell& neighborCell    = *( storage.getCell( edge.neighborCells().at( neighborCellID ) ) );
         auto        cellLocalEdgeID = neighborCell.getLocalEdgeID( edge.getID() );

         const auto basisInCell = algorithms::getMissingIntegersAscending< 2, 4 >(
             {neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
              neighborCell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

         const auto centerIndexInCell = indexing::basisConversion(
             centerIndexOnEdge, basisInCell, {0, 1, 2, 3}, levelinfo::num_microedges_per_edge( level ) );
         const auto cellCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
             edgeCenterOrientation, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );

         for ( const auto& leafOrientationInCell : edgedof::allEdgeDoFOrientations )
         {
            for ( const auto& stencilIt : opr_data[neighborCellID][cellCenterOrientation][leafOrientationInCell] )
            {
               const auto stencilOffset = stencilIt.first;

               const auto leafOrientationOnEdge = edgedof::convertEdgeDoFOrientationCellToFace(
                   leafOrientationInCell, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
               const auto leafIndexInCell = centerIndexInCell + stencilOffset;

               const auto leafIndexOnEdge = indexing::basisConversion(
                   leafIndexInCell, {0, 1, 2, 3}, basisInCell, levelinfo::num_microedges_per_edge( level ) );

               const auto onCellFacesSet = edgedof::macrocell::isOnCellFaces( level, leafIndexInCell, leafOrientationInCell );
               const auto onCellFacesSetOnEdge =
                   edgedof::macrocell::isOnCellFaces( level, leafIndexOnEdge, leafOrientationOnEdge );

               WALBERLA_ASSERT_EQUAL( onCellFacesSet.size(), onCellFacesSetOnEdge.size() )

               uint_t leafArrayIndexOnEdge = std::numeric_limits< uint_t >::max();

               const auto cellLocalIDsOfNeighborFaces =
                   indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID );
               std::vector< uint_t > cellLocalIDsOfNeighborFacesWithLeafOnThem;
               std::set_intersection( cellLocalIDsOfNeighborFaces.begin(),
                                      cellLocalIDsOfNeighborFaces.end(),
                                      onCellFacesSet.begin(),
                                      onCellFacesSet.end(),
                                      std::back_inserter( cellLocalIDsOfNeighborFacesWithLeafOnThem ) );

               if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 0 )
               {
                  // leaf in macro-cell
                  leafArrayIndexOnEdge = edgedof::macroedge::indexOnNeighborCell(
                      level, leafIndexOnEdge.x(), neighborCellID, edge.getNumNeighborFaces(), leafOrientationOnEdge );
               }
               else if ( cellLocalIDsOfNeighborFacesWithLeafOnThem.size() == 1 )
               {
                  // leaf on macro-face
                  WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) )

                  const auto cellLocalFaceID = *cellLocalIDsOfNeighborFacesWithLeafOnThem.begin();
                  const auto facePrimitiveID = neighborCell.neighborFaces().at( cellLocalFaceID );
                  WALBERLA_ASSERT( std::find( edge.neighborFaces().begin(), edge.neighborFaces().end(), facePrimitiveID ) !=
                                   edge.neighborFaces().end() )

                  // The leaf orientation on the edge must be X, Y or XY since it is located on a neighboring face.
                  // Therefore we need to know the three spanning vertex IDs and convert the leaf orientation again.
                  const auto spanningCellLocalVertices = indexing::cellLocalFaceIDsToSpanningVertexIDs.at( cellLocalFaceID );
                  std::array< uint_t, 4 > faceBasisInCell;
                  if ( spanningCellLocalVertices.count( basisInCell.at( 2 ) ) == 1 )
                  {
                     faceBasisInCell = basisInCell;
                  }
                  else
                  {
                     WALBERLA_ASSERT( spanningCellLocalVertices.count( basisInCell.at( 3 ) ) == 1 );
                     faceBasisInCell    = basisInCell;
                     faceBasisInCell[2] = basisInCell.at( 3 );
                     faceBasisInCell[3] = basisInCell.at( 2 );
                  }

                  const auto leafIndexOnEdgeGhostLayer = indexing::basisConversion(
                      leafIndexInCell, {0, 1, 2, 3}, faceBasisInCell, levelinfo::num_microedges_per_edge( level ) );
                  const auto leafOrientationOnEdgeGhostLayer = edgedof::convertEdgeDoFOrientationCellToFace(
                      leafOrientationInCell, faceBasisInCell.at( 0 ), faceBasisInCell.at( 1 ), faceBasisInCell.at( 2 ) );

                  const auto localFaceIDOnEdge = edge.face_index( facePrimitiveID );
                  leafArrayIndexOnEdge         = edgedof::macroedge::indexOnNeighborFace(
                      level, leafIndexOnEdgeGhostLayer.x(), localFaceIDOnEdge, leafOrientationOnEdgeGhostLayer );
               }
               else
               {
                  // leaf on macro-edge
                  WALBERLA_ASSERT_EQUAL( cellLocalIDsOfNeighborFacesWithLeafOnThem.size(), 2 )
                  WALBERLA_ASSERT( !edgedof::macrocell::isInnerEdgeDoF( level, leafIndexInCell, leafOrientationInCell ) )
                  WALBERLA_ASSERT_EQUAL( leafOrientationOnEdge, edgedof::EdgeDoFOrientation::X )
                  leafArrayIndexOnEdge = edgedof::macroedge::index( level, leafIndexOnEdge.x() );
               }
               count++;
            }
         }
      }
   }

   return count;
}

uint_t getEdgeCouplingsTotal( const uint_t&                                                                       level,
                              PrimitiveStorage&                                                                   storage,
                              const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >& operatorID )
{
   std::array< uint_t, 6 > perEdgeCouplings = {0};

   uint_t k = 0;
   for ( PrimitiveID& edgeID : storage.getEdgeIDs() )
   {
      const Edge& edge    = *( storage.getEdge( edgeID ) );
      perEdgeCouplings[k] = getEdgeCouplingsForEdge( level, edge, storage, operatorID );
      k++;
#ifdef SEPARATE_EDGE_COUPLINGS_FOR_EDGES
      WALBERLA_LOG_INFO_ON_ROOT( "couplings for current edge " << k << " = " << perEdgeCouplings[k-1] );
#else
      break;
#endif
   }
   return perEdgeCouplings[0];
}

// ==================================================================================================================================

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   const uint_t maxLevel = 6;

   std::array< uint_t, maxLevel + 1 > nCellTotal;
   std::array< uint_t, maxLevel + 1 > nCell2Cell;
   std::array< uint_t, maxLevel + 1 > nCell2Face;
   std::array< uint_t, maxLevel + 1 > nCell2Edge;

   std::array< uint_t, maxLevel + 1 > nFaceTotal;
   std::array< uint_t, maxLevel + 1 > nFace2Cell;
   std::array< uint_t, maxLevel + 1 > nFace2Face;
   std::array< uint_t, maxLevel + 1 > nFace2Edge;

   std::array< uint_t, maxLevel + 1 > nEdgeTotal;
   std::array< uint_t, maxLevel + 1 > nEdge2Cell;
   std::array< uint_t, maxLevel + 1 > nEdge2Face;
   std::array< uint_t, maxLevel + 1 > nEdge2Edge;

   std::string           meshFileName = "../data/meshes/3D/tet_tilted_1el.msh";
   MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setup( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setup );

   EdgeDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > > massOperEE(
       storage, 0, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "=======================================================" );
   WALBERLA_LOG_INFO_ON_ROOT( "===== Counting Couplings on Single Reference Cell =====" );
   WALBERLA_LOG_INFO_ON_ROOT( "=======================================================" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   // COMPUTATIONS
   for ( uint_t lvl = 0; lvl <= maxLevel; lvl++ )
   {
      uint_t nDoFsOnEdge = numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( lvl );
      uint_t nDoFsOnFace = numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( lvl );

      // get cell couplings
      std::array< uint_t, 3 > count = getCellCouplingsByType( lvl );
      nCell2Cell[lvl]               = count[0];
      nCell2Face[lvl]               = count[1];
      nCell2Edge[lvl]               = count[2];
      nCellTotal[lvl]               = count[0] + count[1] + count[2];

      // get face couplings
      nFace2Cell[lvl] = nCell2Face[lvl] / 4;
      // we disregard couplings to other macro-faces and macro-edges here!
      nFace2Face[lvl] = nDoFsOnFace > 0 ? 5 * nDoFsOnFace - 3 * nDoFsOnEdge : 0;
      nFace2Edge[lvl] = nDoFsOnFace > 0 ? 3 * ( 2 * nDoFsOnEdge - 2 ) : 0;
      // nFaceTotal[lvl] = nFace2Cell[lvl] + nFace2Face[lvl] + nFace2Edge[lvl];
      nFaceTotal[lvl] = getFaceCouplingsTotal( lvl, *storage, massOperEE.getFaceStencil3DID() );

      // get edge couplings
      nEdge2Cell[lvl] = nCell2Edge[lvl] / 6;
      nEdge2Face[lvl] = nDoFsOnFace > 0 ? ( 2 * nDoFsOnEdge - 2 ) : 0;
      nEdge2Edge[lvl] = nDoFsOnEdge; // we disregard couplings to other macro-edges/faces here!
      // nEdgeTotal[lvl] = nEdge2Cell[lvl] + nEdge2Face[lvl] + nEdge2Edge[lvl];
      nEdgeTotal[lvl] = getEdgeCouplingsTotal( lvl, *storage, massOperEE.getEdgeStencil3DID() );
   }

   // REPORTING
   WALBERLA_LOG_INFO_ON_ROOT( "+-----------------------------------------------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "|                    Edge -> Edge                     |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-----------------------------------------------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "|                   Cell Couplings                    |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "| Level |  total  | cell-cell | cell-face | cell-edge |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+" );

   for ( uint_t lvl = 0; lvl <= maxLevel; lvl++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "|   " << lvl << "   | " << std::setw( 7 ) << nCellTotal[lvl] << " | " << std::setw( 9 )
                                        << nCell2Cell[lvl] << " | " << std::setw( 9 ) << nCell2Face[lvl] << " | "
                                        << std::setw( 9 ) << nCell2Edge[lvl] << " | " );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+" );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------------------------------------------------------------------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "|                              Face Couplings                             |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+---------+---------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "| Level |  total  | face-cell | face-face | face-edge | omitted | assumed |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+---------+---------+" );

   for ( uint_t lvl = 0; lvl <= maxLevel; lvl++ )
   {
      uint_t assumed = 15 * ( (1ul << lvl) - 1 );
      WALBERLA_LOG_INFO_ON_ROOT( "|   " << lvl << "   | " << std::setw( 7 ) << nFaceTotal[lvl] << " | " << std::setw( 9 )
                                        << nFace2Cell[lvl] << " | " << std::setw( 9 ) << nFace2Face[lvl] << " | "
                                        << std::setw( 9 ) << nFace2Edge[lvl] << " | " << std::setw( 7 )
                                        << nFaceTotal[lvl] - ( nFace2Cell[lvl] + nFace2Face[lvl] + nFace2Edge[lvl] ) << " | "
                                        << std::setw( 7 ) << assumed << " |" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+---------+---------+" );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------------------------------------------------------------------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "|                               Edge Couplings                            |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+---------+---------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "| Level |  total  | edge-cell | edge-face | edge-edge | omitted | assumed |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+---------+---------+" );

   for ( uint_t lvl = 0; lvl <= maxLevel; lvl++ )
   {
      uint_t assumed = lvl == 0 ? 5 : 4 * ((1ul << (lvl-1)) + 1);

      WALBERLA_LOG_INFO_ON_ROOT( "|   " << lvl << "   | " << std::setw( 7 ) << nEdgeTotal[lvl] << " | " << std::setw( 9 )
                                        << nEdge2Cell[lvl] << " | " << std::setw( 9 ) << nEdge2Face[lvl] << " | "
                                        << std::setw( 9 ) << nEdge2Edge[lvl] << " | " << std::setw( 7 )
                                        << nEdgeTotal[lvl] - ( nEdge2Cell[lvl] + nEdge2Face[lvl] + nEdge2Edge[lvl] ) << " | "
                                        << std::setw( 7 ) << assumed << " |" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+-----------+-----------+-----------+---------+---------+" );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Remark:" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Our computation of the face-* and edge-* couplings neglects certain" );
   WALBERLA_LOG_INFO_ON_ROOT( "couplings. Specifically we omit" );
   WALBERLA_LOG_INFO_ON_ROOT( " (a) (edge-edge) couplings to other macro-edges besides ourself" );
   WALBERLA_LOG_INFO_ON_ROOT( " (b) (edge-face) counts couplings to one face only" );
   WALBERLA_LOG_INFO_ON_ROOT( " (c) (face-face) couplings to other macro-faces besides ourself" );
   WALBERLA_LOG_INFO_ON_ROOT( " (d) (face-edge) couplings to macro-edges not neighbouring the face itself" );
   WALBERLA_LOG_INFO_ON_ROOT( "The omitted values experimentally behave as (assumed)" );
   WALBERLA_LOG_INFO_ON_ROOT( "  * face: 15*(2^level-1)" );
   WALBERLA_LOG_INFO_ON_ROOT( "  * edge: 4*(2^(level-1)+1) - delta(0,level)" );
   WALBERLA_LOG_INFO_ON_ROOT( "Taking into account that each edge couples to two faces, the number of" );
   WALBERLA_LOG_INFO_ON_ROOT( "omitted couplings becomes 6 for all levels > 0 (i.e. one coupling to each of" );
   WALBERLA_LOG_INFO_ON_ROOT( "the four edges sharing a vertex with this edge and one to the two faces" );
   WALBERLA_LOG_INFO_ON_ROOT( "connected to these vertices. On level 0 the number is 5." );

   // COMPARISON TO TRUE VALUES
   std::array< uint_t, 5 > nnzEE = {36, 217, 1474, 10788, 82376};
   // std::array< uint_t, 5 > nnzVE = { 24, 122, 740, 5064, 37264 };

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "COMPARISON:" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "We compare the sum: nCellTotal + 4 * nFaceTotal + 6 * nEdgeTotal to the true nnz" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+---------+------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "| Level |  #NNZ   | our sum | difference |" );
   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+---------+------------+" );

   std::array< uint_t, maxLevel + 1 > sumEdgeDoF2EdgeDoF;
   for ( uint_t lvl = 0; lvl <= 4; lvl++ )
   {
      sumEdgeDoF2EdgeDoF[lvl] = nCellTotal[lvl];
      sumEdgeDoF2EdgeDoF[lvl] += 4 * nFaceTotal[lvl];
      sumEdgeDoF2EdgeDoF[lvl] += 6 * nEdgeTotal[lvl];
      WALBERLA_LOG_INFO_ON_ROOT( "|   " << lvl << "   | " << std::setw( 7 ) << nnzEE[lvl] << " | " << std::setw( 7 )
                                        << sumEdgeDoF2EdgeDoF[lvl] << " | " << std::setw( 10 )
                                        << (int) nnzEE[lvl] - (int) sumEdgeDoF2EdgeDoF[lvl] << " | " );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "+-------+---------+---------+------------+" );
}
