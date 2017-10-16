
#include "tinyhhg_core/p1functionspace/P1DataHandling.hpp"

namespace hhg {

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

std::shared_ptr< VertexP1LocalMatrixMemory > VertexP1LocalMatrixMemoryDataHandling::initialize( const Vertex * const vertex) const
{
  auto vertexP1LocalMatrixMemory = std::make_shared< VertexP1LocalMatrixMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexP1LocalMatrixMemory->addlevel( level, vertex->getNumNeighborFaces());
  }
  return vertexP1LocalMatrixMemory;
}

std::shared_ptr< EdgeP1LocalMatrixMemory > EdgeP1LocalMatrixMemoryDataHandling::initialize( const Edge * const ) const
{
  auto edgeP1LocalMatrixMemory = std::make_shared< EdgeP1LocalMatrixMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeP1LocalMatrixMemory->addlevel( level );
  }
  return edgeP1LocalMatrixMemory;
}

std::shared_ptr< FaceP1LocalMatrixMemory > FaceP1LocalMatrixMemoryDataHandling::initialize( const Face * const ) const
{
  auto faceP1LocalMatrixMemory = std::make_shared< FaceP1LocalMatrixMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceP1LocalMatrixMemory->addlevel( level );
  }
  return faceP1LocalMatrixMemory;
}

}
