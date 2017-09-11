
#include "tinyhhg_core/p1functionspace/P1DataHandling.hpp"

namespace hhg {

std::shared_ptr< VertexP1FunctionMemory< real_t > > VertexP1FunctionMemoryDataHandling::initialize( const Vertex * const vertex ) const
{
  auto vertexP1FunctionMemory = std::make_shared< VertexP1FunctionMemory< real_t > >( vertex->getNumNeighborEdges() );
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexP1FunctionMemory->addlevel( level );
  }
  return vertexP1FunctionMemory;
}

std::shared_ptr< EdgeP1FunctionMemory< real_t > > EdgeP1FunctionMemoryDataHandling::initialize( const Edge * const edge ) const
{
  auto edgeP1FunctionMemory = std::make_shared< EdgeP1FunctionMemory< real_t > >( edge->getNumNeighborFaces() );
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeP1FunctionMemory->addlevel( level );
  }
  return edgeP1FunctionMemory;
}

std::shared_ptr< FaceP1FunctionMemory< real_t > > FaceP1FunctionMemoryDataHandling::initialize( const Face * const ) const
{
  auto faceP1FunctionMemory = std::make_shared< FaceP1FunctionMemory< real_t > >( 0 );
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceP1FunctionMemory->addlevel( level );
  }
  return faceP1FunctionMemory;
}

std::shared_ptr< VertexP1StencilMemory > VertexP1StencilMemoryDataHandling::initialize( const Vertex * const vertex) const
{
  auto vertexP1StencilMemory = std::make_shared< VertexP1StencilMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexP1StencilMemory->addlevel( level, vertex->getNumNeighborEdges());
  }
  return vertexP1StencilMemory;
}

std::shared_ptr< EdgeP1StencilMemory > EdgeP1StencilMemoryDataHandling::initialize( const Edge * const ) const
{
  auto edgeP1StencilMemory = std::make_shared< EdgeP1StencilMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeP1StencilMemory->addlevel( level );
  }
  return edgeP1StencilMemory;
}

std::shared_ptr< FaceP1StencilMemory > FaceP1StencilMemoryDataHandling::initialize( const Face * const ) const
{
  auto faceP1StencilMemory = std::make_shared< FaceP1StencilMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceP1StencilMemory->addlevel( level );
  }
  return faceP1StencilMemory;
}

}
