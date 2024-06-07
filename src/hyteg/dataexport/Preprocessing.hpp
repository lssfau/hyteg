/*
* Copyright (c) 2024 Nils Kohl.
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

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"

#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

#include "sqlite/SQLite.h"

using walberla::double_c;
using walberla::int64_c;
using walberla::int64_t;
using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

/// Interpolates the number of neighboring macro-volumes that a DoF is located on as either a DoF or a ghost-DoF.
///
/// For instance, vertex DoFs that belong to a macro-edge will be assigned the number of macro-volumes that share that edge.
///
/// This is useful since the VTK output (and possibly others) write all DoFs and ghost-DoFs that reside on macro-volumes.
/// Clearly, DoFs at macro-volume boundaries are therefore written multiple times. Sorting that issue out (i.e., writing each DoF
/// exactly once) and maintaining the mesh connectivity information would require global enumeration (which is possible) and
/// sophisticated parallel file writing because the DoFs have to be writing in a specific order _globally_ (which is rather
/// complicated).
///
///
/// This function helps by keeping the simple output routines, but enables adding information about the amount of duplication that
/// can be used for post-processing.
///
/// \tparam ScalarFunctionType function type to interpolate into
/// \param u                   function to interpolate into
/// \param level               refinement level
template < typename ScalarFunctionType >
inline void interpolateNumberOfAdjacentMacroVolumes( ScalarFunctionType& u, uint_t level )
{
   using ValueType = typename ScalarFunctionType::valueType;

   const bool isP1Function = std::is_same_v< typename ScalarFunctionType::Tag, P1FunctionTag >;
   const bool isP2Function = std::is_same_v< typename ScalarFunctionType::Tag, P2FunctionTag >;

   WALBERLA_CHECK( isP1Function || isP2Function,
                   "interpolateNumberOfAdjacentMacroVolumes() only supported for scalar P1 and P2 functions." );

   std::shared_ptr< PrimitiveStorage > storage = u.getStorage();

   for ( auto [pid, vertex] : storage->getVertices() )
   {
      ValueType numNeighborVolumes =
          static_cast< ValueType >( storage->hasGlobalCells() ? vertex->getNumNeighborCells() : vertex->getNumNeighborFaces() );

      WALBERLA_CHECK_GREATER_EQUAL( numNeighborVolumes, 1 );

      if constexpr ( isP1Function )
      {
         auto primitiveDataID = u.getVertexDataID();
         vertexdof::macrovertex::interpolate( level, *vertex, primitiveDataID, numNeighborVolumes );
      }
      else
      {
         auto primitiveDataIDVertexDoF = u.getVertexDoFFunction().getVertexDataID();
         vertexdof::macrovertex::interpolate( level, *vertex, primitiveDataIDVertexDoF, numNeighborVolumes );
      }
   }

   for ( auto [pid, edge] : storage->getEdges() )
   {
      ValueType numNeighborVolumes =
          static_cast< ValueType >( storage->hasGlobalCells() ? edge->getNumNeighborCells() : edge->getNumNeighborFaces() );

      WALBERLA_CHECK_GREATER_EQUAL( numNeighborVolumes, 1 );

      if constexpr ( isP1Function )
      {
         auto primitiveDataID = u.getEdgeDataID();
         vertexdof::macroedge::interpolate( level, *edge, primitiveDataID, numNeighborVolumes );
      }
      else
      {
         auto primitiveDataIDVertexDoF = u.getVertexDoFFunction().getEdgeDataID();
         auto primitiveDataIDEdgeDoF   = u.getEdgeDoFFunction().getEdgeDataID();
         vertexdof::macroedge::interpolate( level, *edge, primitiveDataIDVertexDoF, numNeighborVolumes );
         edgedof::macroedge::interpolate( level, *edge, primitiveDataIDEdgeDoF, numNeighborVolumes );
      }
   }

   for ( auto [pid, face] : storage->getFaces() )
   {
      ValueType numNeighborVolumes =
          static_cast< ValueType >( storage->hasGlobalCells() ? face->getNumNeighborCells() : ValueType( 1 ) );

      if constexpr ( isP1Function )
      {
         auto primitiveDataID = u.getFaceDataID();
         vertexdof::macroface::interpolate( level, *face, primitiveDataID, numNeighborVolumes );
      }
      else
      {
         auto primitiveDataIDVertexDoF = u.getVertexDoFFunction().getFaceDataID();
         auto primitiveDataIDEdgeDoF   = u.getEdgeDoFFunction().getFaceDataID();
         vertexdof::macroface::interpolate( level, *face, primitiveDataIDVertexDoF, numNeighborVolumes );
         edgedof::macroface::interpolate( level, *face, primitiveDataIDEdgeDoF, numNeighborVolumes );
      }
   }

   for ( auto [pid, cell] : storage->getCells() )
   {
      if constexpr ( isP1Function )
      {
         auto primitiveDataID = u.getCellDataID();
         vertexdof::macrocell::interpolate( level, *cell, primitiveDataID, ValueType( 1 ) );
      }
      else
      {
         auto primitiveDataIDVertexDoF = u.getVertexDoFFunction().getCellDataID();
         auto primitiveDataIDEdgeDoF   = u.getEdgeDoFFunction().getCellDataID();
         vertexdof::macrocell::interpolate( level, *cell, primitiveDataIDVertexDoF, ValueType( 1 ) );
         edgedof::macrocell::interpolate( level, *cell, primitiveDataIDEdgeDoF, ValueType( 1 ) );
      }
   }
}

} // namespace hyteg