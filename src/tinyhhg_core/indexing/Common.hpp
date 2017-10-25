
#pragma once

#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {
namespace indexing {

/// Wrapper around Point3D for convenient access to logical indices.
class Index : private PointND< int, 3 >
{
public:
  const int & col() const { return x[0]; }
        int & col()       { return x[0]; }

  const int & row() const { return x[1]; }
        int & row()       { return x[1]; }

  const int & dep() const { return x[2]; }
        int & dep()       { return x[2]; }
};

}
}
