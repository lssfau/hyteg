/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/dgfunctionspace_old/DGFunction.hpp"

#include "hyteg/dgfunctionspace_old/DGFunctionMacroEdge.hpp"
#include "hyteg/dgfunctionspace_old/DGFunctionMacroFace.hpp"
#include "hyteg/dgfunctionspace_old/DGFunctionMacroVertex.hpp"

namespace hyteg {

template < typename ValueType >
void DGFunction_old< ValueType >::projectP1( P1Function< ValueType >& src, uint_t level, DoFType flag, UpdateType updateType )
{
   this->startTiming( "projectP1" );

   src.template startCommunication< Edge, Vertex >( level );
   src.template startCommunication< Face, Edge >( level );
   src.template endCommunication< Edge, Vertex >( level );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         dgfunction::macrovertex::projectP1< ValueType >(
             level, vertex, this->getStorage(), src.getVertexDataID(), this->getVertexDataID(), updateType );
      }
   }

   this->template startCommunication< Vertex, Edge >( level );

   src.template endCommunication< Face, Edge >( level );

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

      if ( testFlag( edgeBC, flag ) )
      {
         dgfunction::macroedge::projectP1< ValueType >(
             level, edge, this->getStorage(), src.getEdgeDataID(), this->getEdgeDataID(), updateType );
      }
   }

   this->template endCommunication< Vertex, Edge >( level );

   this->template startCommunication< Edge, Face >( level );

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         dgfunction::macroface::projectP1< ValueType >(
             level, face, this->getStorage(), src.getFaceDataID(), this->getFaceDataID(), updateType );
      }
   }

   this->template endCommunication< Edge, Face >( level );

   this->stopTiming( "projectP1" );
}

template class DGFunction_old< real_t >;
template class DGFunction_old< int32_t >;
template class DGFunction_old< int64_t >;

} // namespace hyteg
