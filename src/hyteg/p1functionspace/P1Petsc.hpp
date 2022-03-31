/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

// =============
//  P1Operators
// =============
inline void saveIdentityOperator( const P1Function< idx_t >&                  dst,
                                  const std::shared_ptr< SparseMatrixProxy >& mat,
                                  size_t                                      level,
                                  DoFType                                     flag )
{
   const auto storage = dst.getStorage();

   for ( auto& it : storage->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::saveIdentityOperator( vertex, dst.getVertexDataID(), mat, level );
      }
   }

   for ( auto& it : storage->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::saveIdentityOperator( level, edge, dst.getEdgeDataID(), mat );
      }
   }

   if ( level >= 2 )
   {
      for ( auto& it : storage->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::saveIdentityOperator( level, face, dst.getFaceDataID(), mat );
         }
      }

      for ( auto& it : storage->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            vertexdof::macrocell::saveIdentityOperator( level, cell, dst.getCellDataID(), mat );
         }
      }
   }
}

} // namespace hyteg
