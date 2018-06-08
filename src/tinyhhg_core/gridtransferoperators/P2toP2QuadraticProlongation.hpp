
#pragma once

#include "tinyhhg_core/p2functionspace/P2Function.hpp"

namespace hhg {

class P2toP2QuadraticProlongation
{
public:

    inline void operator()( const P2Function <real_t> & function, const uint_t & sourceLevel, const DoFType & flag ) const
    {
      const auto vertexDoFFunction = function.getVertexDoFFunction();
      const auto edgeDoFFunction = function.getEdgeDoFFunction();

      edgeDoFFunction->template communicate< Vertex, Edge >( sourceLevel );
      edgeDoFFunction->template communicate< Edge, Face >( sourceLevel );

      vertexDoFFunction->template communicate< Vertex, Edge >( sourceLevel );
      vertexDoFFunction->template communicate< Edge, Face >( sourceLevel );

      for( const auto& it : function.getStorage()->getFaces() )
      {
        const Face& face = *it.second;

        const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
        if( testFlag( faceBC, flag ) )
        {
          P2::macroface::prolongate< real_t >(
          sourceLevel, face, vertexDoFFunction->getFaceDataID(), edgeDoFFunction->getFaceDataID() );
        }
      }

      for( const auto& it : function.getStorage()->getEdges() )
      {
        const Edge& edge = *it.second;

        const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
        if( testFlag( edgeBC, flag ) )
        {
          P2::macroedge::prolongate< real_t >(
          sourceLevel, edge, vertexDoFFunction->getEdgeDataID(), edgeDoFFunction->getEdgeDataID() );
        }
      }

      for( const auto& it : function.getStorage()->getVertices() )
      {
        const Vertex& vertex = *it.second;

        const DoFType vertexBC = function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
        if( testFlag( vertexBC, flag ) )
        {
          P2::macrovertex::prolongate< real_t >(
          sourceLevel, vertex, vertexDoFFunction->getVertexDataID(), edgeDoFFunction->getVertexDataID() );
        }
      }

    }
};

}
