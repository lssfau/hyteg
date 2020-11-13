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

#include "hyteg/FunctionProperties.hpp"
// #include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Vertex.hpp"

namespace hyteg {
namespace indexing {

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
      WALBERLA_ABORT( "Mission a failure!" );
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
      WALBERLA_ABORT( "Mission a failure!" );
   }

   return nCouplings;
}

template <>
uint_t countLocalDoFCouplings< P2FunctionTag, P2FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   uint_t VV = countLocalDoFCouplings< P1FunctionTag, P1FunctionTag >( storage, level );
   uint_t EE = countLocalDoFCouplings< EdgeDoFFunctionTag, EdgeDoFFunctionTag >( storage, level );
   uint_t EV = countLocalDoFCouplings< EdgeDoFFunctionTag, VertexDoFFunctionTag >( storage, level );

   // WALBERLA_LOG_INFO_ON_ROOT( " V -> V = " << VV );
   // WALBERLA_LOG_INFO_ON_ROOT( " E -> E = " << EE );
   // WALBERLA_LOG_INFO_ON_ROOT( " E -> V = " << EV );

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
uint_t countLocalDoFCouplings< P1VectorFunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
  uint_t dim = storage->hasGlobalCells() ? 3 : 2;
  return dim * countLocalDoFCouplings< P1FunctionTag, P1FunctionTag >( storage, level );
}

template <>
uint_t countLocalDoFCouplings< P1FunctionTag, P2VectorFunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
  uint_t dim = storage->hasGlobalCells() ? 3 : 2;
  return dim * countLocalDoFCouplings< P2FunctionTag, P1FunctionTag >( storage, level );
}

template <>
uint_t countLocalDoFCouplings< P2VectorFunctionTag, P1FunctionTag >( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
  uint_t dim = storage->hasGlobalCells() ? 3 : 2;
  return dim * countLocalDoFCouplings< P2FunctionTag, P1FunctionTag >( storage, level );
}
// -------------------------------------------------------------------------------------------------------------------------------------


} // namespace indexing
} // namespace hyteg
