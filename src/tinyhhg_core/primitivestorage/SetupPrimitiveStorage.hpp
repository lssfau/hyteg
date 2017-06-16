
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/SetupEdge.hpp"
#include "tinyhhg_core/primitives/SetupFace.hpp"
#include "tinyhhg_core/primitives/SetupVertex.hpp"

#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace hhg {

class SetupPrimitiveStorage
{
public:

  SetupPrimitiveStorage( const MeshInfo & meshInfo )
  {

  }

private:

  void generatePrimitiveID();

  std::map< PrimitiveID::IDType, SetupVertex* > vertices_;
  std::map< PrimitiveID::IDType, SetupEdge*  >  edges_;
  std::map< PrimitiveID::IDType, SetupFace* >   faces_;

};

} // namespace hhg
