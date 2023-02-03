/*
 * Copyright (c) 2017-2023 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType >
class FaceDoFFunction_old;

void P1IntegrateDG( FaceDoFFunction_old< real_t >& rhs,
                    P1Function< real_t >&          rhsP1,
                    P1Function< real_t >&          dstP1,
                    uint_t                         level,
                    DoFType                        flag )
{
   rhsP1.startCommunication< Edge, Vertex >( level );
   rhsP1.startCommunication< Face, Edge >( level );

   rhs.startCommunication< Face, Edge >( level );
   rhs.endCommunication< Face, Edge >( level );

   rhs.startCommunication< Edge, Vertex >( level );
   rhs.endCommunication< Edge, Vertex >( level );

   rhsP1.endCommunication< Edge, Vertex >( level );

   auto storage = dstP1.getStorage();

   for ( auto& it : storage->getVertices() )
   {
      Vertex& vertex = *it.second;

      if ( testFlag( dstP1.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macrovertex::integrateDG< real_t >(
             vertex, storage, rhs.getVertexDataID(), rhsP1.getVertexDataID(), dstP1.getVertexDataID(), level );
      }
   }

   dstP1.startCommunication< Vertex, Edge >( level );
   rhsP1.endCommunication< Face, Edge >( level );

   for ( auto& it : storage->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( dstP1.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroedge::integrateDG< real_t >(
             level, edge, storage, rhs.getEdgeDataID(), rhsP1.getEdgeDataID(), dstP1.getEdgeDataID() );
      }
   }

   dstP1.endCommunication< Vertex, Edge >( level );
   dstP1.startCommunication< Edge, Face >( level );

   for ( auto& it : storage->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( dstP1.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         vertexdof::macroface::integrateDG< real_t >(
             level, face, rhs.getFaceDataID(), rhsP1.getFaceDataID(), dstP1.getFaceDataID() );
      }
   }

   dstP1.endCommunication< Edge, Face >( level );
}

} // namespace hyteg
