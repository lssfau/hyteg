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

// Backward Differentiation Formula Scheme
// Assumes that you have a time dependent system of the form
//    dx/dt + Ax = f
// that you want to solve using a BDF of a given order as a first derivative approximation.
// In particular f is assumed to be explicit or a suitable extrapolation ( any implicit parts are assumed to be handled via A ).
// Order 0 implies time independent.
//
// This function expects the time step sizes
// dt0 = t0 - t1, dt1 = t1 - t2, ... to be saved in the stepSizes vector of
// the Function History (state 0 holds dt0, state 1 holds dt1, ...).
// Hence you should call newState( newTimeStepSize ) on your FunctionHistory before trying to solve
// your system and save the new result in state 0.
//
// If not enough past states are available then the function automatically
// applies a BDF scheme of lower order or aborts if the history is empty.
//
// Example: BDF1 (implicit Euler):
//    dx/dt ~ ( x(t^n) - x(t^{n-1}) ) / (t^n - t^{n-1} ) =: x^n/tau - x^{n-1}/tau
// Hence we get:
//    tau * Ax^n + x^n = tau * f^n + x^{n-1}
// Usage:
//    -After you have determined your righthandside f for a given timestep, e.g. by applying some operators and
//     adding up the results use applyRHS to calculate:
//       tau * f + terms depending on past states (i.e. scaled mass applied to past states of x)
//    -At the end of the apply function of your Operator A call applyLHS to calculate:
//       tau * Ax + scaled mass applied to x
template < class DstFunctionType,
           class MassOperatorType,
           class AdditionalDataType = real_t,
           class SrcFunctionType    = DstFunctionType >
class BDFScheme : public TimeDiscretisationScheme< DstFunctionType, MassOperatorType, AdditionalDataType, SrcFunctionType >
{
 public:
   BDFScheme( uint_t order )
   : order_( order )
   {}

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
      uint_t order        = std::min( order_, maximumOrder );

      switch ( order )
      {
      case 0:
         break;
      case 1: {
         const real_t dt0 = history.getStateStepSize( 0 );

         const real_t scaling0 = real_c( 1 );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 2: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );

         const real_t scaling0 = real_c( 1 ) + dt0 / ( dt0 + dt1 );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 3: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );

         const real_t scaling0 = real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 4: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt3 + dt0 + dt1 + dt2 ) );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 5: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );

         const real_t scaling0 = real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) +
                                                 dt0 / ( dt0 + dt1 + dt2 + dt3 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 6: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 7: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 8: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );
         const real_t dt7 = history.getStateStepSize( 7 );

         const real_t scaling0 =
             real_c( 1 ) +
             ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      case 9: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );
         const real_t dt7 = history.getStateStepSize( 7 );
         const real_t dt8 = history.getStateStepSize( 8 );

         const real_t scaling0 =
             real_c( 1 ) +
             ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );

         hyteg::applyGEMV( massOperator, scaling0, src, dt0, dst, level, flag );
      }
      break;
      default:
         WALBERLA_ABORT( "Requested BDF order not implemented!" );
         break;
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
      uint_t order        = std::min( order_, maximumOrder );

      switch ( order )
      {
      case 0:
         break;
      case 1: {
         const real_t dt0 = history.getStateStepSize( 0 );

         const real_t scaling1 = real_c( 1 );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), dt0, dst, level, flag );
      }
      break;
      case 2: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );

         const real_t scaling1 = ( dt0 + dt1 ) / dt1;
         const real_t scaling2 = -( dt0 * dt0 ) / ( dt0 * dt1 + ( dt1 * dt1 ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
      }
      break;
      case 3: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );

         const real_t scaling1 = ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt1 + dt2 ) );
         const real_t scaling2 = -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt0 + dt1 ) * dt2 ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) ) / ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling3, history.getState( 3 ), real_c( 1 ), dst, level, flag );
      }
      break;
      case 4: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );

         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt3 + dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt1 + dt2 ) * ( dt3 + dt1 + dt2 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt3 + dt0 + dt1 + dt2 ) ) / ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt3 + dt2 ) ) );
         const real_t scaling3 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt3 + dt0 + dt1 + dt2 ) ) / ( dt3 * dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) ) /
                                    ( dt3 * ( dt3 + dt2 ) * ( dt3 + dt1 + dt2 ) * ( dt3 + dt0 + dt1 + dt2 ) ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling3, history.getState( 3 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling4, history.getState( 4 ), real_c( 1 ), dst, level, flag );
      }
      break;
      case 5: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );

         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
                                    ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 ) );
         const real_t scaling5 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) ) /
             ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling3, history.getState( 3 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling4, history.getState( 4 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling5, history.getState( 5 ), real_c( 1 ), dst, level, flag );
      }
      break;
      case 6: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );

         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) );
         const real_t scaling4 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
                ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) ) );
         const real_t scaling5 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 );
         const real_t scaling6 = -(
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) ) /
             ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling3, history.getState( 3 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling4, history.getState( 4 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling5, history.getState( 5 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling6, history.getState( 6 ), real_c( 1 ), dst, level, flag );
      }
      break;
      case 7: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );

         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) *
                                   ( dt3 + dt4 + dt5 + dt6 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                    ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 *
                                      ( dt4 + dt5 ) * ( dt4 + dt5 + dt6 ) ) );
         const real_t scaling5 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                 ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 * ( dt5 + dt6 ) );
         const real_t scaling6 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                                    ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                                      ( dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 ) );
         const real_t scaling7 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) ) /
             ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling3, history.getState( 3 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling4, history.getState( 4 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling5, history.getState( 5 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling6, history.getState( 6 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling7, history.getState( 7 ), real_c( 1 ), dst, level, flag );
      }
      break;
      case 8: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );
         const real_t dt7 = history.getStateStepSize( 7 );

         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) );
         const real_t scaling3 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                                 ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) *
                                   ( dt3 + dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) );
         const real_t scaling4 = -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                                      ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                                    ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 *
                                      ( dt4 + dt5 ) * ( dt4 + dt5 + dt6 ) * ( dt4 + dt5 + dt6 + dt7 ) ) );
         const real_t scaling5 = ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                                 ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) *
                                   ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 * ( dt5 + dt6 ) * ( dt5 + dt6 + dt7 ) );
         const real_t scaling6 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
                ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 * ( dt6 + dt7 ) ) );
         const real_t scaling7 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
             ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * dt7 );
         const real_t scaling8 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) ) /
                ( dt7 * ( dt6 + dt7 ) * ( dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling3, history.getState( 3 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling4, history.getState( 4 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling5, history.getState( 5 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling6, history.getState( 6 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling7, history.getState( 7 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling8, history.getState( 8 ), real_c( 1 ), dst, level, flag );
      }
      break;
      case 9: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );
         const real_t dt7 = history.getStateStepSize( 7 );
         const real_t dt8 = history.getStateStepSize( 8 );

         const real_t scaling1 =
             ( ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt1 * ( dt1 + dt2 ) * ( dt1 + dt2 + dt3 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling2 =
             -( ( dt0 * dt0 * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt1 * ( dt0 + dt1 ) * dt2 * ( dt2 + dt3 ) * ( dt2 + dt3 + dt4 ) * ( dt2 + dt3 + dt4 + dt5 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) );
         const real_t scaling3 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt2 * ( dt1 + dt2 ) * ( dt0 + dt1 + dt2 ) * dt3 * ( dt3 + dt4 ) * ( dt3 + dt4 + dt5 ) * ( dt3 + dt4 + dt5 + dt6 ) *
               ( dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling4 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt3 * ( dt2 + dt3 ) * ( dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 ) * dt4 * ( dt4 + dt5 ) *
                  ( dt4 + dt5 + dt6 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 + dt8 ) ) );
         const real_t scaling5 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt4 * ( dt3 + dt4 ) * ( dt2 + dt3 + dt4 ) * ( dt1 + dt2 + dt3 + dt4 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) * dt5 *
               ( dt5 + dt6 ) * ( dt5 + dt6 + dt7 ) * ( dt5 + dt6 + dt7 + dt8 ) );
         const real_t scaling6 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt5 * ( dt4 + dt5 ) * ( dt3 + dt4 + dt5 ) * ( dt2 + dt3 + dt4 + dt5 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * dt6 * ( dt6 + dt7 ) * ( dt6 + dt7 + dt8 ) ) );
         const real_t scaling7 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
             ( dt6 * ( dt5 + dt6 ) * ( dt4 + dt5 + dt6 ) * ( dt3 + dt4 + dt5 + dt6 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) * dt7 * ( dt7 + dt8 ) );
         const real_t scaling8 =
             -( ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) ) /
                ( dt7 * ( dt6 + dt7 ) * ( dt5 + dt6 + dt7 ) * ( dt4 + dt5 + dt6 + dt7 ) * ( dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) *
                  ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) * dt8 ) );
         const real_t scaling9 =
             ( dt0 * dt0 * ( dt0 + dt1 ) * ( dt0 + dt1 + dt2 ) * ( dt0 + dt1 + dt2 + dt3 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) *
               ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) ) /
             ( dt8 * ( dt7 + dt8 ) * ( dt6 + dt7 + dt8 ) * ( dt5 + dt6 + dt7 + dt8 ) * ( dt4 + dt5 + dt6 + dt7 + dt8 ) *
               ( dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) * ( dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) *
               ( dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) * ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );

         dst.assign( { dt0 }, { dst }, level, flag );

         hyteg::applyGEMV( massOperator, scaling1, history.getState( 1 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling2, history.getState( 2 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling3, history.getState( 3 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling4, history.getState( 4 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling5, history.getState( 5 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling6, history.getState( 6 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling7, history.getState( 7 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling8, history.getState( 8 ), real_c( 1 ), dst, level, flag );
         hyteg::applyGEMV( massOperator, scaling9, history.getState( 9 ), real_c( 1 ), dst, level, flag );
      }
      break;
      default:
         WALBERLA_ABORT( "Requested BDF order not implemented!" );
         break;
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
      uint_t order        = std::min( order_, maximumOrder );

      switch ( order )
      {
      case 0:
         break;
      case 1: {
         const real_t dt0 = history.getStateStepSize( 0 );

         const real_t scaling0 = real_c( 1 );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 2: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );

         const real_t scaling0 = real_c( 1 ) + dt0 / ( dt0 + dt1 );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 3: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );

         const real_t scaling0 = real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 4: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt3 + dt0 + dt1 + dt2 ) );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 5: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );

         const real_t scaling0 = real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) +
                                                 dt0 / ( dt0 + dt1 + dt2 + dt3 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 6: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 7: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );

         const real_t scaling0 =
             real_c( 1 ) + ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
                             dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 8: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );
         const real_t dt7 = history.getStateStepSize( 7 );

         const real_t scaling0 =
             real_c( 1 ) +
             ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      case 9: {
         const real_t dt0 = history.getStateStepSize( 0 );
         const real_t dt1 = history.getStateStepSize( 1 );
         const real_t dt2 = history.getStateStepSize( 2 );
         const real_t dt3 = history.getStateStepSize( 3 );
         const real_t dt4 = history.getStateStepSize( 4 );
         const real_t dt5 = history.getStateStepSize( 5 );
         const real_t dt6 = history.getStateStepSize( 6 );
         const real_t dt7 = history.getStateStepSize( 7 );
         const real_t dt8 = history.getStateStepSize( 8 );

         const real_t scaling0 =
             real_c( 1 ) +
             ( dt0 / ( dt0 + dt1 ) + dt0 / ( dt0 + dt1 + dt2 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 ) + dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 ) +
               dt0 / ( dt0 + dt1 + dt2 + dt3 + dt4 + dt5 + dt6 + dt7 + dt8 ) );

         mat->createFromMatrixLinComb( { dt0 }, { mat } );
         hyteg::applyToMatrixScaled( massOperator, scaling0, mat, srcIdx, dstIdx, level, flag );
      }
      break;
      default:
         WALBERLA_ABORT( "Requested BDF order not implemented!" );
         break;
      }
   }

 private:
   uint_t order_; // order of the BFD Scheme
};

} // namespace hyteg