
#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/gridtransferoperators/ProlongationOperator.hpp"

namespace hhg {

class P1toP1QuadraticProlongation : public ProlongationOperator< P1Function< real_t > >
{
public:

    inline void prolongate( const P1Function< real_t > & function, const uint_t & sourceLevel, const DoFType & flag ) const override
    {
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