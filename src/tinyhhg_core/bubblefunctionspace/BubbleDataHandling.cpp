#include "BubbleDataHandling.hpp"

namespace hhg {

VertexBubbleFunctionMemory * VertexBubbleFunctionMemoryDataHandling::initialize( const Vertex * const vertex ) const
{
  VertexBubbleFunctionMemory * vertexBubbleFunctionMemory = new VertexBubbleFunctionMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexBubbleFunctionMemory->addlevel( level, vertex->getNumNeighborEdges() );
  }
  return vertexBubbleFunctionMemory;
}

EdgeBubbleFunctionMemory * EdgeBubbleFunctionMemoryDataHandling::initialize( const Edge * const edge ) const
{
  EdgeBubbleFunctionMemory * edgeBubbleFunctionMemory = new EdgeBubbleFunctionMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeBubbleFunctionMemory->addlevel( level, edge->getNumNeighborFaces() );
  }
  return edgeBubbleFunctionMemory;
}

FaceBubbleFunctionMemory * FaceBubbleFunctionMemoryDataHandling::initialize( const Face * const ) const
{
  FaceBubbleFunctionMemory * faceBubbleFunctionMemory = new FaceBubbleFunctionMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceBubbleFunctionMemory->addlevel( level );
  }
  return faceBubbleFunctionMemory;
}

FaceBubbleStencilMemory * FaceBubbleStencilMemoryDataHandling::initialize( const Face * const ) const
{
  FaceBubbleStencilMemory * faceBubbleStencilMemory = new FaceBubbleStencilMemory();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceBubbleStencilMemory->addlevel( level );
  }
  return faceBubbleStencilMemory;
}

}
