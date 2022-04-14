/*
 * Copyright (c) 2021-2022 Andreas Wagner, Marcus Mohr.
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

#include "core/DataTypes.h"

#include "hyteg/types/types.hpp"

namespace hyteg {

using walberla::real_t;

template < typename Function >
class WeightedJacobiSmoothable
{
 public:
   virtual ~WeightedJacobiSmoothable() = default;

   virtual void smooth_jac( const Function& dst,
                            const Function& rhs,
                            const Function& tmp,
                            real_t          relax,
                            size_t          level,
                            DoFType         flag ) const = 0;
};

template < typename Function >
class ConstantJacobiSmoothable
{
 public:
   virtual ~ConstantJacobiSmoothable() = default;

   virtual void smooth_jac( const Function& dst,
                            const Function& rhs,
                            const Function& tmp,
                            size_t          level,
                            DoFType         flag ) const = 0;
};

template < typename Function >
class GSSmoothable
{
 public:
   virtual ~GSSmoothable() = default;

   virtual void smooth_gs( const Function& dst, const Function& rhs, size_t level, DoFType flag ) const = 0;
};

template < typename Function >
class GSBackwardsSmoothable
{
 public:
   virtual ~GSBackwardsSmoothable() = default;

   virtual void smooth_gs_backwards( const Function& dst, const Function& rhs, size_t level, DoFType flag ) const = 0;
};

template < typename Function >
class SORSmoothable
{
 public:
   virtual ~SORSmoothable() = default;

   virtual void smooth_sor( const Function& dst, const Function& rhs, real_t relax, size_t level, DoFType flag ) const = 0;
};

template < typename Function >
class SORBackwardsSmoothable
{
 public:
   virtual ~SORBackwardsSmoothable() = default;

   virtual void
       smooth_sor_backwards( const Function& dst, const Function& rhs, real_t relax, size_t level, DoFType flag ) const = 0;
};

template < typename Function >
class OperatorWithInverseDiagonal
{
 public:
   virtual ~OperatorWithInverseDiagonal() = default;

   virtual std::shared_ptr< Function > getInverseDiagonalValues() const = 0;
   virtual void computeInverseDiagonalOperatorValues() = 0;
};

} // namespace hyteg
