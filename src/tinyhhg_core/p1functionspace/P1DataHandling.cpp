
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

VertexP1StencilMemory * VertexP1StencilMemoryDataHandling::initialize( const Vertex * const vertex) const
{
  VertexP1StencilMemory * vertexP1StencilMemory = new VertexP1StencilMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexP1StencilMemory->addlevel( level, vertex->getNumNeighborEdges());
  }
  return vertexP1StencilMemory;
}

EdgeP1StencilMemory * EdgeP1StencilMemoryDataHandling::initialize( const Edge * const ) const
{
  EdgeP1StencilMemory * edgeP1StencilMemory = new EdgeP1StencilMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeP1StencilMemory->addlevel( level );
  }
  return edgeP1StencilMemory;
}

FaceP1StencilMemory * FaceP1StencilMemoryDataHandling::initialize( const Face * const ) const
{
  FaceP1StencilMemory * faceP1StencilMemory = new FaceP1StencilMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceP1StencilMemory->addlevel( level );
  }
  return faceP1StencilMemory;
}

}
