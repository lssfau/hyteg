/*
 * Copyright (c) 2020 Daniel Drzisga
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

#include "hyteg/Operator.hpp"
#include "hyteg/composites//P1StokesFunction.hpp"

namespace hyteg {

using walberla::real_t;

class P1ProjectNormalOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   P1ProjectNormalOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const std::function<void(const Point3D&, Point3D& )>& normal_function );

   ~P1ProjectNormalOperator() override = default;

   void apply( const P1StokesFunction< real_t >& dst,
               size_t                            level,
               DoFType                           flag) const;

 private:
   const std::function<void(const Point3D&, Point3D& )> normal_function_;
};

} // namespace hyteg
