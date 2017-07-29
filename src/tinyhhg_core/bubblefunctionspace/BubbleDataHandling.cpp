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

//VertexBubbleStencilMemory * VertexBubbleStencilMemoryDataHandling::initialize( const Vertex * const vertex) const
//{
//  VertexBubbleStencilMemory * vertexBubbleStencilMemory = new VertexBubbleStencilMemory();
//  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
//  {
//    vertexBubbleStencilMemory->addlevel( level, vertex->getNumNeighborEdges());
//  }
//  return vertexBubbleStencilMemory;
//}
//
//EdgeBubbleStencilMemory * EdgeBubbleStencilMemoryDataHandling::initialize( const Edge * const ) const
//{
//  EdgeBubbleStencilMemory * edgeBubbleStencilMemory = new EdgeBubbleStencilMemory();
//  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
//  {
//    edgeBubbleStencilMemory->addlevel( level );
//  }
//  return edgeBubbleStencilMemory;
//}

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
