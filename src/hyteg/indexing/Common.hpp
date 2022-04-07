/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "hyteg/types/pointnd.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {
namespace indexing {

using walberla::uint_t;

class IndexIncrement : public PointND< int, 3 >
{
 public:
   IndexIncrement()
   : PointND< int, 3 >()
   {}
   IndexIncrement( const IndexIncrement& other )
   : PointND< int, 3 >( other )
   {}

   IndexIncrement( const int& x, const int& y, const int& z )
   {
      x_[0] = x;
      x_[1] = y;
      x_[2] = z;
   }

   const int& x() const { return x_[0]; }
   int&       x() { return x_[0]; }

   const int& y() const { return x_[1]; }
   int&       y() { return x_[1]; }

   const int& z() const { return x_[2]; }
   int&       z() { return x_[2]; }

   void setxyz( const int& x, const int& y, const int& z )
   {
      x_[0] = x;
      x_[1] = y;
      x_[2] = z;
   }

   IndexIncrement& operator+=( const IndexIncrement& increment )
   {
      x() += increment.x();
      y() += increment.y();
      z() += increment.z();
      return *this;
   }
};

/// Wrapper around Point3D for convenient access to logical indices.
class Index : public PointND< idx_t, 3 >
{
 public:
   Index()
   : PointND< idx_t, 3 >()
   {}
   Index( const Index& other )
   : PointND< idx_t, 3 >( other )
   {}

   static Index max()
   {
      return Index( std::numeric_limits< idx_t >::max(), std::numeric_limits< idx_t >::max(), std::numeric_limits< idx_t >::max() );
   }

   Index( const idx_t& x, const idx_t& y, const idx_t& z )
   {
      x_[0] = x;
      x_[1] = y;
      x_[2] = z;
   }

   const idx_t& x() const { return x_[0]; }
   idx_t&       x() { return x_[0]; }

   const idx_t& y() const { return x_[1]; }
   idx_t&       y() { return x_[1]; }

   const idx_t& z() const { return x_[2]; }
   idx_t&       z() { return x_[2]; }

   const idx_t& col() const { return x_[0]; }
   idx_t&       col() { return x_[0]; }

   const idx_t& row() const { return x_[1]; }
   idx_t&       row() { return x_[1]; }

   const idx_t& dep() const { return x_[2]; }
   idx_t&       dep() { return x_[2]; }

   Index& operator+=( const IndexIncrement& increment )
   {
      WALBERLA_ASSERT_GREATER_EQUAL( (idx_t) x() + increment.x(), 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( (idx_t) y() + increment.y(), 0 );
      WALBERLA_ASSERT_GREATER_EQUAL( (idx_t) z() + increment.z(), 0 );
      x() += increment.x();
      y() += increment.y();
      z() += increment.z();
      return *this;
   }
};

inline bool operator<( const Index& lhs, const Index& rhs )
{
   return lhs.x() < rhs.x() || ( lhs.x() == rhs.x() && lhs.y() < rhs.y() ) ||
          ( lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() < rhs.z() );
}

inline bool operator<( const IndexIncrement& lhs, const IndexIncrement& rhs )
{
   return lhs.x() < rhs.x() || ( lhs.x() == rhs.x() && lhs.y() < rhs.y() ) ||
          ( lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() < rhs.z() );
}

inline Index operator+( Index lhs, const IndexIncrement& rhs )
{
   lhs += rhs;
   return lhs;
}

inline Index operator+( const IndexIncrement& lhs, Index rhs )
{
   rhs += lhs;
   return rhs;
}

inline IndexIncrement operator+( IndexIncrement lhs, const IndexIncrement& rhs )
{
   lhs += rhs;
   return lhs;
}

inline Index operator*( Index lhs, const idx_t& scalar )
{
   lhs.x() *= scalar;
   lhs.y() *= scalar;
   lhs.z() *= scalar;
   return lhs;
}

inline Index operator*( const idx_t& scalar, Index rhs )
{
   rhs.x() *= scalar;
   rhs.y() *= scalar;
   rhs.z() *= scalar;
   return rhs;
}

inline IndexIncrement operator-( const Index& lhs, const Index& rhs )
{
   return IndexIncrement( (int) lhs.x() - (int) rhs.x(), (int) lhs.y() - (int) rhs.y(), (int) lhs.z() - (int) rhs.z() );
}

inline IndexIncrement operator-( const IndexIncrement& lhs, const IndexIncrement& rhs )
{
   return IndexIncrement( (int) lhs.x() - (int) rhs.x(), (int) lhs.y() - (int) rhs.y(), (int) lhs.z() - (int) rhs.z() );
}

inline bool operator==( const Index& lhs, const Index& rhs )
{
   return lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() == rhs.z();
}

inline bool operator==( const IndexIncrement& lhs, const IndexIncrement& rhs )
{
   return lhs.x() == rhs.x() && lhs.y() == rhs.y() && lhs.z() == rhs.z();
}

inline bool operator!=( const IndexIncrement& lhs, const IndexIncrement& rhs )
{
   return !( lhs == rhs );
}

inline std::ostream& operator<<( std::ostream& os, const Index& index )
{
   os << "( " << index.x() << ", " << index.y() << ", " << index.z() << " )";
   return os;
}

inline std::ostream& operator<<( std::ostream& os, const IndexIncrement& indexIncrement )
{
   os << "( " << indexIncrement.x() << ", " << indexIncrement.y() << ", " << indexIncrement.z() << " )";
   return os;
}

} // namespace indexing
} // namespace hyteg
