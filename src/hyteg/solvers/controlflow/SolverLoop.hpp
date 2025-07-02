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

#include <functional>

#include "core/DataTypes.h"

#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

using walberla::uint_t;

template < typename OperatorType >
class SolverLoop : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   SolverLoop( const std::shared_ptr< Solver< OperatorType > >& solver, const uint_t& iterations )
   : solver_( solver )
   , iterations_( iterations )
   , stopIterationCallback_( []( const OperatorType&, const FunctionType&, const FunctionType&, const uint_t ) { return false; } )
   {}

   SolverLoop( const std::shared_ptr< Solver< OperatorType > >& solver,
               const uint_t&                                    iterations,
               const std::function< bool( const OperatorType&, const FunctionType&, const FunctionType&, const uint_t ) >&                   stopIterationCallback )
   : solver_( solver )
   , iterations_( iterations )
   , stopIterationCallback_( stopIterationCallback )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      stopIterationCallback_( A, x, b, level );
      for ( uint_t i = 0; i < iterations_; i++ )
      {
         solver_->solve( A, x, b, level );
         if ( stopIterationCallback_( A, x, b, level ) )
         {
            break;
         }
      }
   }

   void setStopIterationCallback(
       const std::function< bool( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) >&
           stopIterationCallback )
   {
      stopIterationCallback_ = stopIterationCallback;
   }

 private:
   std::shared_ptr< Solver< OperatorType > > solver_;
   uint_t                                    iterations_;
   std::function< bool( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) >
       stopIterationCallback_;
};

} // namespace hyteg