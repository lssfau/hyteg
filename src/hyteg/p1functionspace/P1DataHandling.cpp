
#include "hyteg/p1functionspace/P1DataHandling.hpp"

namespace hyteg {


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

std::shared_ptr< FaceP1PolynomialMemory > FaceP1PolynomialMemoryDataHandling::initialize( const Face * const ) const
{
  return std::make_shared< FaceP1PolynomialMemory >( );
}

}
