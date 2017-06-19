
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include <map>
#include <vector>

namespace hhg {

using walberla::uint_t;


PrimitiveID PrimitiveStorage::addVertex()
{
  PrimitiveID id(0);
  vertices_[ id.getID() ] = new Vertex( 0, Point3D() );
  return id;
}


} // namespace hhg

