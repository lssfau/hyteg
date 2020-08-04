/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/Algorithms.hpp"
#include "hyteg/LevelWiseMemory.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/Operator.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"

namespace hyteg {

class EdgeDoFProjectNormalOperator final : public Operator< hyteg::EdgeDoFFunction< real_t >, hyteg::EdgeDoFFunction< real_t > >
{
 public:
   EdgeDoFProjectNormalOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                                 size_t                                                   minLevel,
                                 size_t                                                   maxLevel,
                                 const std::function< void( const Point3D&, Point3D& ) >& normal_function );
   ~EdgeDoFProjectNormalOperator() = default;

   void apply( const EdgeDoFFunction< real_t >& dst_u,
               const EdgeDoFFunction< real_t >& dst_v,
               const EdgeDoFFunction< real_t >& dst_w,
               uint_t                           level,
               DoFType                          flag) const;

 private:
   const std::function< void( const Point3D&, Point3D& ) > normal_function_;
};

} // namespace hyteg
