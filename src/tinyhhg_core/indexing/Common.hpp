
#pragma once

#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {
namespace indexing {

using walberla::uint_t;

class IndexIncrement : protected PointND< int, 3 >
{
public:

  IndexIncrement()                               : PointND< int, 3 >()        {}
  IndexIncrement( const IndexIncrement & other ) : PointND< int, 3 >( other ) {}

  IndexIncrement( const int & x, const int & y, const int & z )
  {
    x_[0] = x;
    x_[1] = y;
    x_[2] = z;
  }

  const int & x() const { return x_[0]; }
        int & x()       { return x_[0]; }

  const int & y() const { return x_[1]; }
        int & y()       { return x_[1]; }

  const int & z() const { return x_[2]; }
        int & z()       { return x_[2]; }
};


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

  Index & operator+=( const IndexIncrement & increment )
  {
    WALBERLA_ASSERT_GREATER_EQUAL( (int)x() + increment.x(), 0 );
    WALBERLA_ASSERT_GREATER_EQUAL( (int)y() + increment.y(), 0 );
    WALBERLA_ASSERT_GREATER_EQUAL( (int)z() + increment.z(), 0 );
    x() += increment.x();
    y() += increment.y();
    z() += increment.z();
    return *this;
  }

};


inline Index operator+( Index lhs, const IndexIncrement & rhs )
{
  lhs += rhs;
  return lhs;
}

inline Index operator+( const IndexIncrement & lhs, Index rhs )
{
  rhs += lhs;
  return rhs;
}

}
}
