/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/EmptySolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

namespace hyteg {

template < class OperatorType, class pressureBlockPreconditionerType >
class StokesBlockDiagonalPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   StokesBlockDiagonalPreconditioner(
       const std::shared_ptr< PrimitiveStorage >&                                    storage,
       uint_t                                                                        minLevel,
       uint_t                                                                        maxLevel,
       uint_t                                                                        velocityPreconditionSteps,
       std::shared_ptr< hyteg::Solver< typename OperatorType::VelocityOperator_T > > velocityBlockPreconditioner =
           std::make_shared< hyteg::IdentityPreconditioner< typename OperatorType::VelocityOperator_T > >() )
   : velocityPreconditionSteps_( velocityPreconditionSteps )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary )
   , velocityBlockPreconditioner_( velocityBlockPreconditioner )
   , pressureBlockPreconditioner_( std::make_shared< pressureBlockPreconditionerType >( storage, minLevel, maxLevel ) )
   {}

   // y = M^{-1} * x
   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, uint_t level ) override
   {
      for ( uint_t steps = 0; steps < velocityPreconditionSteps_; steps++ )
      {
         for ( uint_t k = 0; k < x.uvw().getDimension(); k++ )
         {
            velocityBlockPreconditioner_->solve( A.getA(), x.uvw()[k], b.uvw()[k], level );
         }

         pressureBlockPreconditioner_->apply( b.p(), x.p(), level, flag_, Replace );
      }
   }

 private:
   uint_t                                                                        velocityPreconditionSteps_;
   hyteg::DoFType                                                                flag_;
   std::shared_ptr< hyteg::Solver< typename OperatorType::VelocityOperator_T > > velocityBlockPreconditioner_;
   std::shared_ptr< pressureBlockPreconditionerType >                            pressureBlockPreconditioner_;
};

} // namespace hyteg
