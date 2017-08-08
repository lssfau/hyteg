
#include <tinyhhg_core/primitives/Primitive.hpp>
#include <tinyhhg_core/primitivestorage/PrimitiveStorage.hpp>

namespace hhg {

void Primitive::getNeighborPrimitives( std::vector< PrimitiveID > & neighborPrimitives ) const
{
  getNeighborVertices( neighborPrimitives );

  std::vector< PrimitiveID > neighborEdges;
  std::vector< PrimitiveID > neighborFaces;

  getNeighborEdges( neighborEdges );
  getNeighborFaces( neighborFaces );

  neighborPrimitives.insert( neighborPrimitives.end(), neighborEdges.begin(), neighborEdges.end() );
  neighborPrimitives.insert( neighborPrimitives.end(), neighborFaces.begin(), neighborFaces.end() );

}

} // namespace hhg

