
#pragma once

#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {
namespace indexing {

using walberla::uint_t;

/// Wrapper around Point3D for convenient access to logical indices.
class Index : protected PointND< uint_t, 3 >
{
public:

  Index()                      : PointND< uint_t, 3 >()        {}
  Index( const Index & other ) : PointND< uint_t, 3 >( other ) {}

  Index( const uint_t & x, const uint_t & y, const uint_t & z )
  {
    x_[0] = x;
    x_[1] = y;
    x_[2] = z;
  }

  const uint_t & x() const { return x_[0]; }
        uint_t & x()       { return x_[0]; }

  const uint_t & y() const { return x_[1]; }
        uint_t & y()       { return x_[1]; }

  const uint_t & z() const { return x_[2]; }
        uint_t & z()       { return x_[2]; }

  const uint_t & col() const { return x_[0]; }
        uint_t & col()       { return x_[0]; }

  const uint_t & row() const { return x_[1]; }
        uint_t & row()       { return x_[1]; }

  const uint_t & dep() const { return x_[2]; }
        uint_t & dep()       { return x_[2]; }
};

}
}
