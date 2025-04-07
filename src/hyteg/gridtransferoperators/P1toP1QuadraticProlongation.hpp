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

#pragma once

#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {

class P1toP1QuadraticProlongation : public ProlongationOperator< P1Function< real_t > >
{
public:

    inline void prolongate( const P1Function< real_t > & function, const uint_t & sourceLevel, const DoFType & flag ) const override
    {
      WALBERLA_CHECK( !function.getStorage()->hasGlobalCells(), "P1 quadratic prolongation not implemented in 3D." )

      const uint_t destinationLevel = sourceLevel + 1;

      for ( const auto & it : function.getStorage()->getVertices() )
      {
        const Vertex & vertex = *it.second;

        if ( testFlag( function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
        {
          const auto srcData = vertex.getData( function.getVertexDataID())->getPointer( sourceLevel );
          auto dstData = vertex.getData( function.getVertexDataID())->getPointer( destinationLevel );
          prolongateMacroVertex( srcData, dstData, sourceLevel );
        }
      }

      function.startCommunication<Vertex, Edge>( destinationLevel );

      for ( const auto & it : function.getStorage()->getEdges() )
      {
        const Edge & edge = *it.second;

        if ( testFlag( function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
        {
          const auto srcData = edge.getData( function.getEdgeDataID())->getPointer( sourceLevel );
          auto dstData = edge.getData( function.getEdgeDataID())->getPointer( destinationLevel );
          prolongateMacroEdge( srcData, dstData, sourceLevel );
        }
      }

      function.endCommunication<Vertex, Edge>( destinationLevel );
      function.startCommunication<Edge, Face>( destinationLevel );

      for ( const auto& it : function.getStorage()->getFaces() )
      {
        const Face & face = *it.second;

        if ( testFlag( function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
        {
          const auto srcData = face.getData( function.getFaceDataID())->getPointer( sourceLevel );
          auto dstData = face.getData( function.getFaceDataID())->getPointer( destinationLevel );
          prolongateMacroFace( srcData, dstData, sourceLevel );
        }
      }

      function.endCommunication<Edge, Face>( destinationLevel );
    }

private:

    void prolongateMacroVertex( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const;

    void prolongateMacroEdge( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const;

    void prolongateMacroFace( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const;

};

}
