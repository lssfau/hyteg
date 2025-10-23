/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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
#include "core/config/Config.h"

#include "hyteg/functions/FunctionHistory.hpp"
#include "hyteg/operators/GEMV.hpp"
#include "hyteg/types/types.hpp"

#include "TimeDiscretisationScheme.hpp"

namespace hyteg {

// Crank-Nicolson time discretisation
// Assumes that you have a time dependent system of the form
//    dx/dt + Ax = f
// that you want to solve using the Crank-Nicolson method.
// In particular f is assumed to be explicit or a suitable extrapolation ( any implicit parts are assumed to be handled via A ).
//
// This function expects the time step size
// dt0 = t0 - t1 to be saved in the stepSizes vector of
// the Function History (state 0 holds dt0, state 1 holds dt1, ...).
// Hence you should call newState( newTimeStepSize ) on your FunctionHistory before trying to solve
// your system and save the new result in state 0.
//
// If no past state is available then the function calls WALBERLA_ABORT.
//
// The right-hand-side at the last time step is provided by the user via a FEM function fLast.
//
// Scheme:
//    Let tau := tn - t{n-1} be the time step.
//    x^n + tau * 1/2 * Ax^n  = tau * 1/2 * ( f + fLast ) + x^{n-1}
// Note: Following the scheme above fLast also needs to contain -Ax^{n-1}!
//
// Usage:
//    -After you have determined your righthandside f for a given timestep, e.g. by applying some operators and
//     adding up the results and you have determined fLast (you probably have to assign f after the last solve)
//     use applyRHS to calculate:
//       tau * 1/2 * ( f + fLast ) + term depending on the past state ( i.e. mass applied to x^{n-1} )
//    -At the end of the apply function of your Operator A call applyLHS to calculate:
//       tau * 1/2 * Ax^n + term depending on the current state ( i.e. mass applied to x^n )
template < class DstFunctionType,
           class MassOperatorType,
           class AdditionalDataType = real_t,
           class SrcFunctionType    = DstFunctionType >
class CrankNicolsonScheme
: public TimeDiscretisationScheme< DstFunctionType, MassOperatorType, AdditionalDataType, SrcFunctionType >
{
 public:
   CrankNicolsonScheme( const std::shared_ptr< DstFunctionType >& fLast )
   : fLast_( fLast )
   {}

   void                               setFLast( const std::shared_ptr< DstFunctionType >& fLast ) { fLast_ = fLast; }
   std::shared_ptr< DstFunctionType > getFLast() { return fLast_; }

   // Warning! Only apply this AFTER you have completely determined the rest of the LHS!
   void applyLHS( const FunctionHistory< SrcFunctionType, AdditionalDataType >& history,
                  const SrcFunctionType&                                        src,
                  const DstFunctionType&                                        dst,
                  const MassOperatorType&                                       massOperator,
                  const uint_t                                                  level,
                  const hyteg::DoFType                                          flag ) override
   {
      uint_t numberOfInitialisedStates = history.getNumberOfInitialisedStates();
      if ( numberOfInitialisedStates == 0 )
      {
         WALBERLA_ABORT( "Empty function history!" );
      }
      uint_t maximumOrder = numberOfInitialisedStates - 1;

      // technically you don't need the last time step here, however we still throw an error
      // to alert the user about potentially improper usage
      if ( maximumOrder >= 1 )
      {
         const real_t dt0 = history.getStateStepSize( 0 );

         const real_t scaling0 = real_c( 1 );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0 * real_c( 0.5 ), dst, level, flag );
      }
      else
      {
         WALBERLA_ABORT( "No past state in function history!" );
      }
   }

   // Warning! Only apply this AFTER you have completely determined the rest of the RHS!
   void applyRHS( const FunctionHistory< SrcFunctionType, AdditionalDataType >& history,
                  const DstFunctionType&                                        dst,
                  const MassOperatorType&                                       massOperator,
                  const uint_t                                                  level,
                  const hyteg::DoFType                                          flag ) override
   {
      uint_t numberOfInitialisedStates = history.getNumberOfInitialisedStates();
      if ( numberOfInitialisedStates == 0 )
      {
         WALBERLA_ABORT( "Empty function history!" );
      }
      uint_t maximumOrder = numberOfInitialisedStates - 1;

      if ( maximumOrder >= 1 )
      {
         const real_t dt0 = history.getStateStepSize( 0 );

         const real_t scaling1 = real_c( 1 );

         // calculate the rhs average
         dst.assign( { real_c( 0.5 ), real_c( 0.5 ) }, { dst, *fLast_ }, level, hyteg::All );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), dt0, dst, level, flag );
      }
      else
      {
         WALBERLA_ABORT( "No past state in function history!" );
      }
   }

   // Warning! Only apply this AFTER you have completely assembled the rest of the matrix!
   void toMatrix( const FunctionHistory< SrcFunctionType, AdditionalDataType >&          history,
                  const MassOperatorType&                                                massOperator,
                  const std::shared_ptr< hyteg::SparseMatrixProxy >&                     mat,
                  const typename SrcFunctionType::template FunctionType< hyteg::idx_t >& srcIdx,
                  const typename DstFunctionType::template FunctionType< hyteg::idx_t >& dstIdx,
                  uint_t                                                                 level,
                  hyteg::DoFType                                                         flag ) override
   {
      uint_t numberOfInitialisedStates = history.getNumberOfInitialisedStates();
      if ( numberOfInitialisedStates == 0 )
      {
         WALBERLA_ABORT( "Empty function history!" );
      }
      uint_t maximumOrder = numberOfInitialisedStates - 1;

      // technically you don't need the last time step here, however we still throw an error
      // to alert the user about potentially improper usage
      if ( maximumOrder >= 1 )
      {
         const real_t dt0 = history.getStateStepSize( 0 );

         const real_t scaling0 = real_c( 1 );

         mat->createFromMatrixLinComb( { dt0 * real_c( 0.5 ) }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      else
      {
         WALBERLA_ABORT( "No past state in function history!" );
      }
   }

 private:
   std::shared_ptr< DstFunctionType > fLast_;
};

} // namespace hyteg