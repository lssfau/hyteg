
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

  IndexIncrement & operator+=( const IndexIncrement & increment )
  {
    x() += increment.x();
    y() += increment.y();
    z() += increment.z();
    return *this;
  }

  void serialize( walberla::mpi::SendBuffer & sendBuffer ) const
  {
    for ( size_t index = 0; index < 2; index++ )
    {
      sendBuffer << x_[index];
    }
  }

  void deserialize( walberla::mpi::RecvBuffer & recvBuffer )
  {
    for ( size_t index = 0; index < 2; index++ )
    {
      recvBuffer >> x_[index];
    }
  }

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

inline bool operator< ( const Index & lhs, const Index & rhs )
{
  return lhs.x() < rhs.x() || ( lhs.x() == rhs.x() && lhs.y() < rhs.y() ) || ( lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() < rhs.z() );
}

inline bool operator< ( const IndexIncrement & lhs, const IndexIncrement & rhs )
{
  return lhs.x() < rhs.x() || ( lhs.x() == rhs.x() && lhs.y() < rhs.y() ) || ( lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() < rhs.z() );
}

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

inline IndexIncrement operator+( IndexIncrement lhs, const IndexIncrement & rhs )
{
  lhs += rhs;
  return lhs;
}

inline Index operator*( Index lhs, const uint_t & scalar )
{
  lhs.x() *= scalar;
  lhs.y() *= scalar;
  lhs.z() *= scalar;
  return lhs;
}

inline Index operator*( const uint_t & scalar, Index rhs )
{
  rhs.x() *= scalar;
  rhs.y() *= scalar;
  rhs.z() *= scalar;
  return rhs;
}

inline IndexIncrement operator-( const Index & lhs, const Index & rhs )
{
  return IndexIncrement( (int) lhs.x() - (int) rhs.x(), (int) lhs.y() - (int) rhs.y(), (int) lhs.z() - (int) rhs.z() );
}

inline IndexIncrement operator-( const IndexIncrement & lhs, const IndexIncrement & rhs )
{
  return IndexIncrement( (int) lhs.x() - (int) rhs.x(), (int) lhs.y() - (int) rhs.y(), (int) lhs.z() - (int) rhs.z() );
}

inline bool operator==( const Index & lhs, const Index & rhs )
{
  return lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() == rhs.z();
}

inline bool operator==( const IndexIncrement & lhs, const IndexIncrement & rhs )
{
  return lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() == rhs.z();
}

inline std::ostream & operator<<( std::ostream & os, const Index & index )
{
  os << "( " << index.x() << ", " << index.y() << ", " << index.z() << " )";
  return os;
}

inline std::ostream & operator<<( std::ostream & os, const IndexIncrement & indexIncrement )
{
  os << "( " << indexIncrement.x() << ", " << indexIncrement.y() << ", " << indexIncrement.z() << " )";
  return os;
}


/// Mapping from [-1, 0, 1]^3 to [0, ..., 26]
inline uint_t mapDirectionToArrayIndex( const IndexIncrement & indexIncrement )
{
  WALBERLA_ASSERT_LESS_EQUAL( indexIncrement.x(), 1 );
  WALBERLA_ASSERT_LESS_EQUAL( indexIncrement.y(), 1 );
  WALBERLA_ASSERT_LESS_EQUAL( indexIncrement.z(), 1 );

  WALBERLA_ASSERT_GREATER_EQUAL( indexIncrement.x(), -1 );
  WALBERLA_ASSERT_GREATER_EQUAL( indexIncrement.y(), -1 );
  WALBERLA_ASSERT_GREATER_EQUAL( indexIncrement.z(), -1 );

  const uint_t bitMask = (1 << (indexIncrement.x() + 7)) | (1 << (indexIncrement.y() + 4)) | (1 << (indexIncrement.z() + 1));
  switch ( bitMask )
  {
  case 0b001001001: // (-1, -1, -1)
    return 0;
  case 0b001001010: // (-1, -1, 0)
    return 1;
  case 0b001001100: // (-1, -1, 1)
    return 2;
  case 0b001010001: // (-1, 0, -1)
    return 3;
  case 0b001010010: // (-1, 0, 0)
    return 4;
  case 0b001010100: // (-1, 0, 1)
    return 5;
  case 0b001100001: // (-1, 1, -1)
    return 6;
  case 0b001100010: // (-1, 1, 0)
    return 7;
  case 0b001100100: // (-1, 1, 1)
    return 8;
  case 0b010001001: // (0, -1, -1)
    return 9;
  case 0b010001010: // (0, -1, 0)
    return 10;
  case 0b010001100: // (0, -1, 1)
    return 11;
  case 0b010010001: // (0, 0, -1)
    return 12;
  case 0b010010010: // (0, 0, 0)
    return 13;
  case 0b010010100: // (0, 0, 1)
    return 14;
  case 0b010100001: // (0, 1, -1)
    return 15;
  case 0b010100010: // (0, 1, 0)
    return 16;
  case 0b010100100: // (0, 1, 1)
    return 17;
  case 0b100001001: // (1, -1, -1)
    return 18;
  case 0b100001010: // (1, -1, 0)
    return 19;
  case 0b100001100: // (1, -1, 1)
    return 20;
  case 0b100010001: // (1, 0, -1)
    return 21;
  case 0b100010010: // (1, 0, 0)
    return 22;
  case 0b100010100: // (1, 0, 1)
    return 23;
  case 0b100100001: // (1, 1, -1)
    return 24;
  case 0b100100010: // (1, 1, 0)
    return 25;
  case 0b100100100: // (1, 1, 1)
    return 26;
  default:
    WALBERLA_ASSERT( false, "Invalid stencil direction" );
    return std::numeric_limits< uint_t >::max();
  }
}


}
}

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
typename G    // Growth policy of SendBuffer
>
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const hhg::indexing::IndexIncrement & indexIncrement )
{
  indexIncrement.serialize( buf );
  return buf;
}

template< typename T // Element type  of RecvBuffer
>
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, hhg::indexing::IndexIncrement & indexIncrement )
{
  indexIncrement.deserialize( buf );
  return buf;
}

}
}