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
#pragma once

#include "hyteg/facedofspace/FaceDoFFunction.hpp"
#include "hyteg/dgfunctionspace/DGFunctionMacroEdge.hpp"
#include "hyteg/dgfunctionspace/DGFunctionMacroFace.hpp"
#include "hyteg/dgfunctionspace/DGFunctionMacroVertex.hpp"

namespace hyteg {

template < typename ValueType >
class DGFunction : public FaceDoFFunction< ValueType >
{
 public:

   DGFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : FaceDoFFunction< ValueType >( name, storage, minLevel, maxLevel )
   {}

   inline void projectP1( P1Function< real_t >& src, uint_t level, DoFType flag, UpdateType updateType = Replace );
};

template < typename ValueType >
void DGFunction< ValueType >::projectP1( P1Function< real_t >& src, uint_t level, DoFType flag, UpdateType updateType )
{
   this->startTiming( "projectP1" );

   src.startCommunication< Edge, Vertex >( level );
   src.startCommunication< Face, Edge >( level );
   src.endCommunication< Edge, Vertex >( level );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         dgfunction::macrovertex::projectP1< real_t >(
             level, vertex, this->getStorage(), src.getVertexDataID(), this->getVertexDataID(), updateType );
      }
   }

   this->template startCommunication< Vertex, Edge >( level );

   src.endCommunication< Face, Edge >( level );

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

      if ( testFlag( edgeBC, flag ) )
      {
         dgfunction::macroedge::projectP1< real_t >(
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
         dgfunction::macroface::projectP1< real_t >(
             level, face, this->getStorage(), src.getFaceDataID(), this->getFaceDataID(), updateType );
      }
   }

   this->template endCommunication< Edge, Face >( level );

   this->stopTiming( "projectP1" );
}

} // namespace hyteg
