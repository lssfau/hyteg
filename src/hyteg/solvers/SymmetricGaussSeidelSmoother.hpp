/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl
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

#include "hyteg/solvers/Solver.hpp"
#include "hyteg/types/flags.hpp"

namespace hyteg {

template < class OperatorType >
class SymmetricGaussSeidelSmoother : public Solver< OperatorType >
{
 public:
   SymmetricGaussSeidelSmoother()
   : flag_( hyteg::Inner | hyteg::NeumannBoundary )
   {}

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      A.smooth_gs( x, b, level, flag_ );
      A.smooth_gs_backwards( x, b, level, flag_ );
   }

 private:
   DoFType flag_;
};

} // namespace hyteg