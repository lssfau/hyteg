
#pragma once

#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {
namespace indexing {

using walberla::uint_t;

/// Wrapper around Point3D for convenient access to logical indices.
class Index : private PointND< uint_t, 3 >
{
public:
  const uint_t & col() const { return x[0]; }
        uint_t & col()       { return x[0]; }

  const uint_t & row() const { return x[1]; }
        uint_t & row()       { return x[1]; }

  const uint_t & dep() const { return x[2]; }
        uint_t & dep()       { return x[2]; }
};

}
}
