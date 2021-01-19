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

#include "hyteg/indexing/CouplingCount.hpp"

#include "hyteg/functions/FunctionProperties.hpp"
// #include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Vertex.hpp"

namespace hyteg {
namespace indexing {

/// Get information on EdgeDoF couplings of a macro-cell on certain refinement level
///
/// The function returns the number of couplings from interior EdgeDoFs of a macro-cell to
///
/// 0: other interior EdgeDoFs
/// 1: EdgeDoFs on the four attached macro-faces of the cell
/// 2: EdgeDoFs on the six attached macro-edges of the cell
///
/// for the specified refinement level. Note that the numbers are hardcoded up to level 7.
/// If you need different levels run the refCellCouplingCount app and add the values.
std::array< uint_t, 3 > getEdgeDoFToEdgeDoFCouplingsForCell( uint_t level )
{
   std::array< uint_t, 8 > cell2cell = {0, 1, 298, 5244, 58304, 543176, 4676568, 38787064};
   std::array< uint_t, 8 > cell2face = {0, 12, 228, 1524, 7572, 33492, 140628, 576084};
   std::array< uint_t, 8 > cell2edge = {0, 0, 12, 36, 84, 180, 372, 756};

   if ( level <= 7 )
   {
      return {cell2cell[level], cell2face[level], cell2edge[level]};
   }
   else
   {
      WALBERLA_ABORT( "getEdgeDoFToEdgeDoFCouplingsForCell: values for requested level not available!" );
   }
}

/// Get information on EdgeDoF to VertexDoF couplings of a macro-cell on certain refinement level
///
/// The function returns the sum of the couplings from interior VertexDoFs of a macro-cell to
///
/// 0: interior EdgeDoFs of the cell
/// 1: EdgeDoFs on the four attached macro-faces of the cell
/// 2: EdgeDoFs on the six attached macro-edges of the cell (no. is always 0)
///
/// for the specified refinement level.
///
/// Note that the numbers are hardcoded up to level 7.
/// If you need different levels run the refCellCouplingCount app and add the values.
uint_t getEdgeDoFToVertexDoFCouplingsForCell( uint_t level )
{
   std::array< uint_t, 8 > cell2cell = {0, 0, 30, 1450, 20930, 216050, 1947730, 16511250};
   std::array< uint_t, 8 > cell2face = {0, 0, 20, 300, 1820, 8700, 37820, 157500};

   if ( level <= 7 )
   {
      return cell2cell[level] + cell2face[level];
   }
   else
   {
      WALBERLA_ABORT( "getEdgeDoFToVertexDoFCouplingsForCell: values for requested level not available!" );
   }
}

/// Get information on EdgeDoF to VertexDoF couplings of a macro-face on certain refinement level
///
/// The function returns the number of the couplings from interior VertexDoFs of a macro-face to
/// interior EdgeDoFs of the cell for the specified refinement level.
///
/// Note that the numbers are hardcoded up to level 7.
/// If you need different levels run the refCellCouplingCount app and add the values.
uint_t getEdgeDoFToVertexDoFCouplingsForFaceToCell( uint_t level )
{
   std::array< uint_t, 8 > face2cell = {0, 0, 37, 339, 1855, 8535, 36487, 150759};

   if ( level <= 7 )
   {
      return face2cell[level];
   }
   else
   {
      WALBERLA_ABORT( "getEdgeDoFToVertexDoFCouplingsForFaceToCell: values for requested level not available!" );
   }
}

/// Get information on EdgeDoF to VertexDoF couplings of a macro-edge on certain refinement level
///
/// The function returns the number of the couplings from interior VertexDoFs of a macro-edge to
/// interior EdgeDoFs of one attached cell for the specified refinement level.
///
/// Note that the numbers are hardcoded up to level 7.
/// If you need different levels run the refCellCouplingCount app and add the values.
uint_t getEdgeDoFToVertexDoFCouplingsForEdgeToCell( uint_t level, uint_t logicalEdgeID )
{
   std::array< std::array< uint_t, 6 >, 7 > edge2cell;
   edge2cell[0] = {0, 0, 0, 0, 0, 0};
   edge2cell[1] = {1, 1, 1, 1, 1, 1};
   edge2cell[2] = {7, 15, 7, 7, 15, 7};
   edge2cell[3] = {19, 43, 19, 19, 43, 19};
   edge2cell[4] = {43, 99, 43, 43, 99, 43};
   edge2cell[5] = {91, 211, 91, 91, 211, 91};
   edge2cell[6] = {187, 435, 187, 187, 435, 187};
   edge2cell[7] = {379, 883, 379, 379, 883, 379};
   if ( level <= 7 )
   {
      return edge2cell[level][logicalEdgeID];
   }
   else
   {
      WALBERLA_ABORT( "getEdgeDoFToVertexDoFCouplingsForFaceToCell: values for requested level not available!" );
   }
}

std::vector< PrimitiveID > getEdgesAttachedToEdge( const std::shared_ptr< PrimitiveStorage >& storage, const PrimitiveID& edgeID )
{
   const Edge* edge = storage->getEdge( edgeID );

   std::vector< PrimitiveID > nbrEdges;

   std::set< PrimitiveID > vEdges;
   std::set< PrimitiveID > cEdges;

   // find all edges attached to the two vertices
   // that are neighbour vertices of the edge
   for ( auto vtxID : edge->neighborVertices() )
   {
      const Vertex* vtx = storage->getVertex( vtxID );
      for ( auto nbrEdgeID : vtx->neighborEdges() )
      {
         vEdges.insert( nbrEdgeID );
      }
   }
   vEdges.erase( edgeID );

   // find all edges attached to the cell neighbours
   // of the edge
   for ( auto cellID : edge->neighborCells() )
   {
      const Cell* cell = storage->getCell( cellID );
      for ( auto nbrEdgeID : cell->neighborEdges() )
      {
         cEdges.insert( nbrEdgeID );
      }
   }
   cEdges.erase( edgeID );

   // compute the intersection (without edge itself)
   for ( auto element : vEdges )
   {
      if ( cEdges.count( element ) == 1 )
         nbrEdges.push_back( element );
   }

   return nbrEdges;
}

template <>
uint_t countLocalDoFCouplings< P1FunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   uint_t nCouplings = 0;

   // ============
   //  2D SETTING
   // ============
   if ( !storage->hasGlobalCells() )
   {
      // for P1-P1 we have 7 couplings for each inner face and edge dof
      nCouplings += 7 * numberOfInnerDoFs< P1FunctionTag, Face >( level ) * storage->getFaceIDs().size();

      // now we add the edge couplings
      std::vector< PrimitiveID > edgeIDs = storage->getEdgeIDs();
      for ( auto edgeID : edgeIDs )
      {
         const Edge* edge         = storage->getEdge( edgeID );
         uint_t      nConnections = 2 * edge->getNumNeighborFaces() + 3;
         nCouplings += nConnections * numberOfInnerDoFs< P1FunctionTag, Edge >( level );
      }

      // now we add the vertex couplings
      std::vector< PrimitiveID > vertexIDs = storage->getVertexIDs();
      for ( auto vtxID : vertexIDs )
      {
         const Vertex* vtx = storage->getVertex( vtxID );
         nCouplings += 1 + vtx->getNumNeighborEdges();
      }
   }

   // ============
   //  3D SETTING
   // ============
   else
   {
      // for P1-P1 we have 15 couplings for each inner cell dof
      nCouplings += 15 * numberOfInnerDoFs< P1FunctionTag, Cell >( level ) * storage->getCellIDs().size();

      // now we add the face couplings
      std::vector< PrimitiveID > faceIDs = storage->getFaceIDs();
      for ( auto faceID : faceIDs )
      {
         const Face* face         = storage->getFace( faceID );
         uint_t      nConnections = 4 * face->getNumNeighborCells() + 7;
         nCouplings += nConnections * numberOfInnerDoFs< P1FunctionTag, Face >( level );
      }

      // now we add the edge couplings
      if ( level >= 1 )
      {
         std::vector< PrimitiveID > edgeIDs = storage->getEdgeIDs();
         for ( auto edgeID : edgeIDs )
         {
            const Edge* edge = storage->getEdge( edgeID );

            // couplings on faces
            uint_t nConnections = 2 * edge->getNumNeighborFaces() + 3;
            nCouplings += nConnections * numberOfInnerDoFs< P1FunctionTag, Edge >( level );

            // couplings into cell interiors:
            //
            // we have one of these, if the edge has index 2 or 4 as seen from the perspective
            // of the neighbouring tetrahedron
            for ( const auto& macroCellID : edge->neighborCells() )
            {
               const auto   macroCell   = storage->getCell( macroCellID );
               const uint_t localEdgeID = macroCell->getLocalEdgeID( edgeID );
               if ( localEdgeID == 2 || localEdgeID == 4 )
               {
                  nCouplings += numberOfInnerDoFs< P1FunctionTag, Edge >( level );
               }
            }
         }
      }

      // Now add vertex couplings
      std::vector< PrimitiveID > vertexIDs = storage->getVertexIDs();
      for ( auto vtxID : vertexIDs )
      {
         const Vertex* vtx = storage->getVertex( vtxID );
         nCouplings += 1 + vtx->getNumNeighborEdges();
      }
   }

   return nCouplings;
}

template <>
uint_t countLocalDoFCouplings< EdgeDoFFunctionTag, EdgeDoFFunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                         uint_t                                     level )
{
   uint_t nCouplings = 0;

   // ============
   //  2D SETTING
   // ============
   if ( !storage->hasGlobalCells() )
   {
      // add macro-face couplings
      nCouplings += 5 * numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level ) * storage->getFaceIDs().size();

      // add macro-edge couplings
      std::vector< PrimitiveID > edgeIDs = storage->getEdgeIDs();
      for ( auto edgeID : edgeIDs )
      {
         const Edge* edge         = storage->getEdge( edgeID );
         uint_t      nConnections = 2 * edge->getNumNeighborFaces() + 1;
         nCouplings += nConnections * numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level );
      }
   }

   // ============
   //  3D SETTING
   // ============
   else
   {
      // -----------------------------------------------------------------------
      //  determine basic coupling numbers for a single cell as building blocks
      // -----------------------------------------------------------------------

      uint_t nDoFsOnEdge = numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level );
      uint_t nDoFsOnFace = numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level );

      // CELL
      std::array< uint_t, 3 > cellCouplings = getEdgeDoFToEdgeDoFCouplingsForCell( level );
      uint_t                  nCell2Cell    = cellCouplings[0];
      uint_t                  nCell2Face    = cellCouplings[1];
      uint_t                  nCell2Edge    = cellCouplings[2];

      // FACE
      // we include couplings to other macro-faces and macro-edges not attached
      // to this macro-face into nFace2Cell (to avoid double counting)
      uint_t nFace2Cell = nCell2Face / 4 + 15 * ( ( 1ul << level ) - 1 );
      uint_t nFace2Face = nDoFsOnFace > 0 ? 5 * nDoFsOnFace - 3 * nDoFsOnEdge : 0;
      uint_t nFace2Edge = nDoFsOnFace > 0 ? 3 * ( 2 * nDoFsOnEdge - 2 ) : 0;

      // EDGE
      // we include the couplings to faces and edges that are only connected to
      // an edge by its vertices into nEdge2Cell (to avoid double counting)
      uint_t nEdge2Cell = nCell2Edge / 6;
      uint_t nEdge2Face = nDoFsOnFace > 0 ? ( 2 * nDoFsOnEdge - 2 ) : 0;
      uint_t nEdge2Edge = nDoFsOnEdge;
      // WALBERLA_LOG_INFO_ON_ROOT( "nEdge2Cell = " << nEdge2Cell );
      // WALBERLA_LOG_INFO_ON_ROOT( "nEdge2Face = " << nEdge2Face );
      // WALBERLA_LOG_INFO_ON_ROOT( "nEdge2Edge = " << nEdge2Edge );

      // --------------------------------------------------
      //  now sum up couplings for the existing primitives
      // --------------------------------------------------

      // CELLS
      uint_t nCells = storage->getCellIDs().size();
      nCouplings += nCells * ( nCell2Cell + nCell2Face + nCell2Edge );
      // WALBERLA_LOG_INFO_ON_ROOT( " --> nCouplings CELL (stop) = " << nCouplings );

      // FACES
      for ( auto faceID : storage->getFaceIDs() )
      {
         const Face* face = storage->getFace( faceID );
         nCouplings += face->getNumNeighborCells() * nFace2Cell + nFace2Face + nFace2Edge;
      }
      // WALBERLA_LOG_INFO_ON_ROOT( " --> nCouplings FACE (stop) = " << nCouplings );

      // EDGES
      for ( auto edgeID : storage->getEdgeIDs() )
      {
         const Edge* edge = storage->getEdge( edgeID );
         nCouplings += nEdge2Edge;
         nCouplings += edge->getNumNeighborFaces() * nEdge2Face;
         nCouplings += edge->getNumNeighborCells() * nEdge2Cell;

         // one extra coupling for each macro-edge that connects to this edge besides ourselves
         nCouplings += getEdgesAttachedToEdge( storage, edgeID ).size();

         if( level == 0 ) {
           nCouplings += 1 * edge->getNumNeighborCells();  // for opposite edge in cell
         }
         else {
           nCouplings += 2 * edge->getNumNeighborCells();  // for two non-attached faces of each attached cell
         }
      }

      // WALBERLA_LOG_INFO_ON_ROOT( " --> nCouplings EDGE (stop) = " << nCouplings );
   }

   return nCouplings;
}

template <>
uint_t countLocalDoFCouplings< EdgeDoFFunctionTag, VertexDoFFunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                           uint_t                                     level )
{
   uint_t nCouplings = 0;

   // ============
   //  2D SETTING
   // ============
   if ( !storage->hasGlobalCells() )
   {
      // add macro-face couplings
      nCouplings += 12 * numberOfInnerDoFs< VertexDoFFunctionTag, Face >( level ) * storage->getFaceIDs().size();

      // add macro-edge couplings
      std::vector< PrimitiveID > edgeIDs = storage->getEdgeIDs();
      for ( auto edgeID : edgeIDs )
      {
         const Edge* edge         = storage->getEdge( edgeID );
         uint_t      nConnections = 5 * edge->getNumNeighborFaces() + 2;
         nCouplings += nConnections * numberOfInnerDoFs< VertexDoFFunctionTag, Edge >( level );
      }

      // add macro-vertex couplings
      std::vector< PrimitiveID > vertexIDs = storage->getVertexIDs();
      for ( auto vtxID : vertexIDs )
      {
         const Vertex* vtx = storage->getVertex( vtxID );
         nCouplings += vtx->getNumNeighborEdges() + vtx->getNumNeighborFaces();
      }
   }

   // ============
   //  3D SETTING
   // ============
   else
   {
      uint_t eDoFsOnEdge = numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level );
      uint_t vDoFsOnEdge = numberOfInnerDoFs< VertexDoFFunctionTag, Edge >( level );
      uint_t vDoFsOnFace = numberOfInnerDoFs< VertexDoFFunctionTag, Face >( level );

      uint_t nCellTotal = getEdgeDoFToVertexDoFCouplingsForCell( level );

      uint_t nFace2Cell = getEdgeDoFToVertexDoFCouplingsForFaceToCell( level );
      uint_t nFace2Edge = level < 2 ? 0 : eDoFsOnEdge - 2;
      uint_t nFace2Face = level < 2 ? 0 : 12 * vDoFsOnFace - 3 * nFace2Edge;
      uint_t nFace2Lost = level < 2 ? 0 : 10 * ( ( 1u << level ) - 2 );

      uint_t nEdge2Edge = 2 * vDoFsOnEdge;
      uint_t nEdge2Face = level == 0 ? 0 : 5 * vDoFsOnEdge - 2;

      // CELLS
      uint_t nCells = storage->getCellIDs().size();
      nCouplings += nCellTotal * nCells;
      nCouplings += 4 * ( nFace2Lost + nFace2Cell ) * nCells;
      // WALBERLA_LOG_INFO_ON_ROOT( " --> nCouplings Cell (stop) = " << nCouplings );

      // FACES
      nCouplings += (nFace2Face + 3 * nFace2Edge) * storage->getFaceIDs().size();

      // WALBERLA_LOG_INFO_ON_ROOT( " --> nCouplings Face (stop) = " << nCouplings );

      // EDGES

      // edges are tricky because some coupling numbers depend on the local index
      // of the edge with respect to the attached cell
      if ( level > 0 )
      {
         for ( auto edgeID : storage->getEdgeIDs() )
         {
            const Edge* edge = storage->getEdge( edgeID );
            nCouplings += nEdge2Edge;
            nCouplings += edge->getNumNeighborFaces() * nEdge2Face;

            // one coupling for each macro-edge that connects to this edge besides ourself
            if ( level > 0 )
              nCouplings += getEdgesAttachedToEdge( storage, edgeID ).size();

            for ( auto cellID : edge->neighborCells() )
            {
               const Cell*  nbrCell     = storage->getCell( cellID );
               const uint_t localEdgeID = nbrCell->getLocalEdgeID( edgeID );
               nCouplings += getEdgeDoFToVertexDoFCouplingsForEdgeToCell( level, localEdgeID );
               if ( localEdgeID == 1 || localEdgeID == 4 )
               {
                  nCouplings += 6;
               }
               else
               {
                  nCouplings += 2;
               }
            }
         }
      }
      // WALBERLA_LOG_INFO_ON_ROOT( " --> nCouplings Edges (stop) = " << nCouplings );

      // VERTICES
      for ( auto vtxID : storage->getVertexIDs() )
      {
         const Vertex*           vtx = storage->getVertex( vtxID );
         std::set< PrimitiveID > edgesInSupport;
         for ( auto cellID : vtx->neighborCells() )
         {
            const Cell* nbrCell = storage->getCell( cellID );
            for ( auto edgeID : nbrCell->neighborEdges() )
            {
               edgesInSupport.insert( edgeID );
            }
         }
         nCouplings += edgesInSupport.size();
      }
      // WALBERLA_LOG_INFO_ON_ROOT( " --> nCouplings Vertex (stop) = " << nCouplings );
   }

   // HACK
   // nCouplings += level == 2 ? 3 : 0;    <-- so the error seems to be with V->E

   return nCouplings;
}

template <>
uint_t countLocalDoFCouplings< P2FunctionTag, P2FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   uint_t VV = countLocalDoFCouplings< P1FunctionTag, P1FunctionTag >( storage, level );
   uint_t EE = countLocalDoFCouplings< EdgeDoFFunctionTag, EdgeDoFFunctionTag >( storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( " V -> V = " << VV );
   WALBERLA_LOG_INFO_ON_ROOT( " E -> E = " << EE );
   uint_t EV = countLocalDoFCouplings< EdgeDoFFunctionTag, VertexDoFFunctionTag >( storage, level );
   WALBERLA_LOG_INFO_ON_ROOT( " E -> V = " << EV );

   return VV + EE + 2 * EV;
}

template <>
uint_t countLocalDoFCouplings< P2FunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   uint_t nCouplings = 0;

   // ============
   //  2D SETTING
   // ============
   if ( !storage->hasGlobalCells() )
   {
      // add macro-face couplings
      nCouplings += 19 * numberOfInnerDoFs< P1FunctionTag, Face >( level ) * storage->getFaceIDs().size();

      // add macro-edge couplings
      std::vector< PrimitiveID > edgeIDs = storage->getEdgeIDs();
      for ( auto edgeID : edgeIDs )
      {
         const Edge* edge         = storage->getEdge( edgeID );
         uint_t      nConnections = 7 * edge->getNumNeighborFaces() + 5;
         nCouplings += nConnections * numberOfInnerDoFs< VertexDoFFunctionTag, Edge >( level );
      }

      // add vertex couplings
      std::vector< PrimitiveID > vertexIDs = storage->getVertexIDs();
      for ( auto vtxID : vertexIDs )
      {
         const Vertex* vtx = storage->getVertex( vtxID );
         nCouplings += 1 + 2 * vtx->getNumNeighborEdges() + 1 * vtx->getNumNeighborFaces();
      }
   }

   // ============
   //  3D SETTING
   // ============
   else
   {
      WALBERLA_ABORT( "Mission a failure!" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "p2 <-> p1 = " << nCouplings );

   return nCouplings;
}

// -------------------------------------------------------------------------------------------------------------------------------------
// NOTE: Should we ever have 2D surface operators in 3D the implementations below might become a problem,
//       as they equates dimension of the vector function with that of the mesh
template <>
uint_t countLocalDoFCouplings< P1VectorFunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                     uint_t                                     level )
{
   uint_t dim = storage->hasGlobalCells() ? 3 : 2;
   return dim * countLocalDoFCouplings< P1FunctionTag, P1FunctionTag >( storage, level );
}

template <>
uint_t countLocalDoFCouplings< P1FunctionTag, P2VectorFunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                     uint_t                                     level )
{
   uint_t dim = storage->hasGlobalCells() ? 3 : 2;
   return dim * countLocalDoFCouplings< P2FunctionTag, P1FunctionTag >( storage, level );
}

template <>
uint_t countLocalDoFCouplings< P2VectorFunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                     uint_t                                     level )
{
   uint_t dim = storage->hasGlobalCells() ? 3 : 2;
   return dim * countLocalDoFCouplings< P2FunctionTag, P1FunctionTag >( storage, level );
}
// -------------------------------------------------------------------------------------------------------------------------------------

} // namespace indexing
} // namespace hyteg
