/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr.
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

#include "hyteg/forms/Form.hpp"

namespace hyteg {

class P2Form : public Form
{
 public:
   virtual ~P2Form() {}

   // 2D P2 VertexDoF
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   // 3D P2 VertexDoF
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
  
  /// Describe a degree of freedom in a P2 element on a tetrahedron
  /// by two vertex indices. There are two cases:
  ///
  /// (a) Both vertex indices are identical, then this is the
  ///     index of the dof associated with this vertex.
  ///
  /// (b) The two vertex indices are different, then this is the
  ///     index of the dof associated with the midpoint of the
  ///     tet's edge given by those two vertices.
  typedef std::array<uint_t,2> dofPosByVertexPair3D;

};

} // namespace hyteg
