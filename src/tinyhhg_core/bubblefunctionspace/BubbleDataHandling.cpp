#include "BubbleDataHandling.hpp"

namespace hhg {

std::shared_ptr< VertexBubbleFunctionMemory > VertexBubbleFunctionMemoryDataHandling::initialize( const Vertex * const vertex ) const
{
  auto vertexBubbleFunctionMemory = std::make_shared< VertexBubbleFunctionMemory >( vertex->getNumNeighborFaces() );
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexBubbleFunctionMemory->addlevel( level );
  }
  return vertexBubbleFunctionMemory;
}

std::shared_ptr< EdgeBubbleFunctionMemory > EdgeBubbleFunctionMemoryDataHandling::initialize( const Edge * const edge ) const
{
  auto edgeBubbleFunctionMemory = std::make_shared< EdgeBubbleFunctionMemory >( edge->getNumNeighborFaces() );
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeBubbleFunctionMemory->addlevel( level );
  }
  return edgeBubbleFunctionMemory;
}

std::shared_ptr< FaceBubbleFunctionMemory > FaceBubbleFunctionMemoryDataHandling::initialize( const Face * const ) const
{
  auto faceBubbleFunctionMemory = std::make_shared< FaceBubbleFunctionMemory >( 0 );
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceBubbleFunctionMemory->addlevel( level );
  }
  return faceBubbleFunctionMemory;
}

std::shared_ptr< FaceBubbleStencilMemory > FaceBubbleStencilMemoryDataHandling::initialize( const Face * const ) const
{
  auto faceBubbleStencilMemory = std::make_shared< FaceBubbleStencilMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceBubbleStencilMemory->addlevel( level );
  }
  return faceBubbleStencilMemory;
}

}
