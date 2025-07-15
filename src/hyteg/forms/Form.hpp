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

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/types/PointND.hpp"
#include "hyteg/types/Matrix.hpp"

namespace hyteg {

/// Base class for all forms
class Form
{
 public:
   Form()
   : geometryMap_( std::make_shared< IdentityMap >() )
   {}

   virtual ~Form() {}

   /// Set the geometry/blending map for the form
   ///
   /// \note
   /// - This method is used e.g. by the ElementwiseOperators.
   /// - In the case of the FEniCS forms the map is ignored.
   virtual void setGeometryMap( const std::shared_ptr< GeometryMap >& geometryMap ) const
   {
      // Check would make sense. However, there are some corner cases, where
      // in 3D nobody calls setGeometryMap on the P2RowSumForm, in which case
      // that will pass a nullptr on to an underlying FenicsForm. Which is
      // sort of okay as the latter will ignore the map anyway :(
      // WALBERLA_ASSERT_NOT_NULLPTR( geometryMap );
      geometryMap_ = geometryMap;
   }

 protected:
   mutable std::shared_ptr< GeometryMap > geometryMap_;
};

} // namespace hyteg
