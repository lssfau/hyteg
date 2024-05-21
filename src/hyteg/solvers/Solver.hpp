/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

namespace hyteg {

template < class OperatorType >
class Solver
{
 public:
   /// solves the system A * x = b
   virtual void solve( const OperatorType&             A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const walberla::uint_t         level ) = 0;

   virtual ~Solver() = default;
};
} // namespace hyteg

// namespace evostencils
namespace evostencils {

enum SolverType 
{
  Empty = 0,
  CG = 1,
  Chebyshev = 2,
  GaussSeidel = 3,
  GKB = 4,
  GMRES= 5,
  Minres = 6,
  SOR = 7,
  StokesPCG = 8,
  SymmetricGaussSeidel = 9,
  SymmetricSOR = 10,
  Uzawa = 11,
  WeightedJacobi = 12,
};
} // namespace evostencils