
#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"

namespace hhg {

class P1toP1LinearRestriction
{
public:

    inline void operator()( const P1Function< real_t > & function, const uint_t & sourceLevel, const DoFType & flag ) const
    {
      const uint_t destinationLevel = sourceLevel - 1;
      const auto storage = function.getStorage();
      const auto boundaryCondition = function.getBoundaryCondition();

      // start pulling vertex halos
      function.startCommunication< Edge, Vertex >( sourceLevel );

      // start pulling edge halos
      function.startCommunication< Face, Edge >( sourceLevel );

      // end pulling vertex halos
      function.endCommunication< Edge, Vertex >( sourceLevel );

      for ( const auto & it : storage->getVertices())
      {
        const Vertex & vertex = *it.second;
        const auto srcData = vertex.getData( function.getVertexDataID())->getPointer( sourceLevel );
        auto dstData = vertex.getData( function.getVertexDataID())->getPointer( destinationLevel );

        if ( testFlag( boundaryCondition.getBoundaryType( vertex.getMeshBoundaryFlag()), flag ))
        {
          restrictMacroVertex( srcData, dstData, sourceLevel, vertex.getNumNeighborEdges() );
        }
      }

      function.startCommunication< Vertex, Edge >( destinationLevel );

      // end pulling edge halos
      function.endCommunication< Face, Edge >( sourceLevel );

      for ( const auto & it : storage->getEdges())
      {
        const Edge & edge = *it.second;
        const auto srcData = edge.getData( function.getEdgeDataID())->getPointer( sourceLevel );
        auto dstData = edge.getData( function.getEdgeDataID())->getPointer( destinationLevel );

        if ( testFlag( boundaryCondition.getBoundaryType( edge.getMeshBoundaryFlag()), flag ))
        {
          restrictMacroEdge( srcData, dstData, sourceLevel, edge.getNumNeighborFaces() );
        }
      }

      function.endCommunication< Vertex, Edge >( destinationLevel );

      function.startCommunication< Edge, Face >( destinationLevel );

      for ( const auto & it : storage->getFaces())
      {
        const Face & face = *it.second;
        const auto srcData = face.getData( function.getFaceDataID())->getPointer( sourceLevel );
        auto dstData = face.getData( function.getFaceDataID())->getPointer( destinationLevel );

        if ( testFlag( boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag()), flag ))
        {
          restrictMacroFace( srcData, dstData, sourceLevel, 0 );
        }
      }

      function.endCommunication< Edge, Face >( destinationLevel );

    }

private:

    void restrictMacroVertex( const real_t *src, real_t *dst, const uint_t & sourceLevel, const uint_t & numNeighborEdges ) const;

    void restrictMacroEdge( const real_t *src, real_t *dst, const uint_t & sourceLevel, const uint_t & numNeighborFaces ) const;

    void restrictMacroFace( const real_t *src, real_t *dst, const uint_t & sourceLevel, const uint_t & numNeighborCells ) const;

};

}