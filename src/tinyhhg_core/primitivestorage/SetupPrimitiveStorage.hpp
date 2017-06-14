
#pragma once

#include "core/debug/Debug.h"
#include "tinyhhg_core/types/pointnd.hpp"

#include <map>

namespace hhg {

using walberla::uint_t;
using walberla::real_t;

class MeshInfo
{
public:

  MeshInfo( const std::string & meshFileName );

private:

  std::map< uint_t, Point3D > vertices_;

};


class SetupPrimitiveStorage
{
public:

  SetupPrimitiveStorage( const MeshInfo & meshInfo )
  {

  }

};

} // namespace hhg
