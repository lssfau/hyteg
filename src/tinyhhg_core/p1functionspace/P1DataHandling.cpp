
#include "tinyhhg_core/p1functionspace/P1DataHandling.hpp"

namespace hhg {

VertexP1FunctionMemory * VertexP1FunctionMemoryDataHandling::initialize( const Vertex * const vertex ) const
{
  VertexP1FunctionMemory * vertexP1FunctionMemory = new VertexP1FunctionMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexP1FunctionMemory->addlevel( level, vertex->getNumNeighborEdges() );
  }
  return vertexP1FunctionMemory;
}

EdgeP1FunctionMemory * EdgeP1FunctionMemoryDataHandling::initialize( const Edge * const edge ) const
{
  EdgeP1FunctionMemory * edgeP1FunctionMemory = new EdgeP1FunctionMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeP1FunctionMemory->addlevel( level, edge->getNumNeighborFaces() );
  }
  return edgeP1FunctionMemory;
}

FaceP1FunctionMemory * FaceP1FunctionMemoryDataHandling::initialize( const Face * const ) const
{
  FaceP1FunctionMemory * faceP1FunctionMemory = new FaceP1FunctionMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceP1FunctionMemory->addlevel( level );
  }
  return faceP1FunctionMemory;
}

}
