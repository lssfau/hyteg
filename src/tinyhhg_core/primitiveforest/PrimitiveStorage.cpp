
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "tinyhhg_core/primitiveforest/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"

#include <map>
#include <vector>

namespace walberla {
namespace hhg {

  PrimitiveID PrimitiveStorage::addPrimitive()
  {
    PrimitiveID id;
    primitives_[ id.getID() ] = new Primitive();
    return id;
  }

  const Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID & id ) const
  {
    return NULL;
  }

  Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID & id )
  {
    return primitives_[ id.getID() ];
  }

  bool PrimitiveStorage::primitiveExistsLocally( const PrimitiveID & id ) const
  {
    return false;
  }


} // namespace hhg
} // namespace walberla
