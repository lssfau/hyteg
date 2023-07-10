/*
 * Copyright (c) 2022 Daniel Bauer.
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
#include "core/timing/TimingTree.h"

#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace hyteg {
namespace n1e1 {

using walberla::real_t;
using walberla::uint_t;

template < class N1E1OperatorType, class P1LaplaceOperatorType >
class HybridSmoother : public Solver< N1E1OperatorType >
{
 public:
   /// \brief Hiptmair's hybrid smoother for curl-curl problems.
   ///
   /// For details see: R. Hiptmair, "Multigrid Method for Maxwell’s Equations," [10.1137/S0036142997326203](https://doi.org/10.1137/S0036142997326203).
   ///
   /// \param storage           A PrimitiveStorage instance.
   /// \param p1LaplaceOperator The operator in potential space.
   /// \param n1e1Smoother      A Solver instance that smoothes error components in the complement of the null space
   ///                          of the curl operator.
   /// \param p1Smoother        A Solver instance that smoothes error components in the null space of the curl operator.
   /// \param minLevel          Minimum level on which memory for temporary functions is allocated.
   /// \param maxLevel          Maximum level on which memory for temporary functions is allocated.
   /// \param n1e1SmoothSteps   The number of smoothing steps in the complement of the null space of the curl operator.
   /// \param p1SmoothSteps     The number of smoothing steps in the null space of the curl operator.
   ///
   HybridSmoother( const std::shared_ptr< PrimitiveStorage >&         storage,
                   std::shared_ptr< P1LaplaceOperatorType >           p1LaplaceOperator,
                   std::shared_ptr< Solver< N1E1OperatorType > >      n1e1Smoother,
                   std::shared_ptr< Solver< P1LaplaceOperatorType > > p1Smoother,
                   uint_t                                             minLevel,
                   uint_t                                             maxLevel,
                   uint_t                                             n1e1SmoothSteps = 1,
                   uint_t                                             p1SmoothSteps   = 1 );

   void solve( const N1E1OperatorType&             A,
               const N1E1VectorFunction< real_t >& x,
               const N1E1VectorFunction< real_t >& b,
               const uint_t                        level ) override;

 private:
   uint_t n1e1SmoothSteps_;
   uint_t p1SmoothSteps_;

   hyteg::DoFType flag_;

   std::shared_ptr< P1LaplaceOperatorType >           p1LaplaceOperator_;
   std::shared_ptr< Solver< N1E1OperatorType > >      n1e1Smoother_;
   std::shared_ptr< Solver< P1LaplaceOperatorType > > p1Smoother_;

   N1E1VectorFunction< real_t > vectorResidual_;
   P1Function< real_t >         scalarResidual_;
   P1Function< real_t >         scalarPotential_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

} // namespace n1e1
} // namespace hyteg
