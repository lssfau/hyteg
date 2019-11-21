/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#include "core/timing/Timer.h"

#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

template < typename OperatorType >
class TimedSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   TimedSolver( const std::shared_ptr< Solver< OperatorType > >& solver )
   : solver_( solver )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timer_.start();
      solver_->solve( A, x, b, level );
      timer_.end();
   }

   walberla::WcTimer& getTimer() { return timer_; }

 private:
   std::shared_ptr< Solver< OperatorType > > solver_;
   walberla::WcTimer                          timer_;
};

} // namespace hyteg