
#pragma once

#include "primitives/face.hpp"
#include "primitives/edge.hpp"

#include <vector>

namespace walberla {
namespace hhg {

class PrimitiveStorage : private walberla::NonCopyable
{
public:

  /// Iterators for traversing all locally allocated faces, edges and vertices
  std::iterator beginFaces()    { return faces_.begin(); }
  std::iterator endFaces()      { return faces_.end(); }
  std::iterator beginEdges()    { return edges_.begin(); }
  std::iterator endEdges()      { return edges_.end(); }
  std::iterator beginVertices() { return edges_.begin(); }
  std::iterator endVertices()   { return edges_.end(); }

private:

  std::map< PrimitiveID::IDType, ::hhg::Face* >   faces_;
  std::map< PrimitiveID::IDType, ::hhg::Edge* >   edges_;
  std::map< PrimitiveID::IDType, ::hhg::Vertex* > vertices_;

};

} // namespace hhg
} // namespace walberla

