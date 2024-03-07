/*
 * Copyright (c) 2023 Andreas Burkhart.
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

#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/solvers/Solver.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

template < class OperatorType, class SubstituteOperatorType >
class SubstitutePreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   SubstitutePreconditioner( std::shared_ptr< Solver< SubstituteOperatorType > > substituteSolver,
                             std::shared_ptr< SubstituteOperatorType >           substituteOperator )
   : holdsReference_( false )
   , substituteSolver_( substituteSolver )
   , substituteOperator_( substituteOperator )
   , substituteOperatorRef_( *substituteOperator )
   {}

   SubstitutePreconditioner( std::shared_ptr< Solver< SubstituteOperatorType > > substituteSolver,
                             SubstituteOperatorType&                             substituteOperator )
   : holdsReference_( true )
   , substituteSolver_( substituteOperator )
   , substituteOperatorRef_( substituteOperator )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_UNUSED( A );
      if ( holdsReference_ )
      {
         substituteSolver_->solve( substituteOperatorRef_.get(), x, b, level );
      }
      else
      {
         substituteSolver_->solve( *substituteOperator_, x, b, level );
      }
   }

 private:
   const bool                                          holdsReference_;
   std::shared_ptr< Solver< SubstituteOperatorType > > substituteSolver_;
   std::shared_ptr< SubstituteOperatorType >           substituteOperator_;
   std::reference_wrapper< SubstituteOperatorType >    substituteOperatorRef_;
};

template < class OperatorType, class SubstituteOperatorType >
class VectorSubstitutePreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   VectorSubstitutePreconditioner( std::shared_ptr< Solver< SubstituteOperatorType > > substituteSolver,
                                   std::shared_ptr< SubstituteOperatorType >           substituteOperator )
   : holdsReference_( false )
   , substituteSolver_( substituteSolver )
   , substituteOperator_( substituteOperator )
   , substituteOperatorRef_( *substituteOperator )
   {}

   VectorSubstitutePreconditioner( std::shared_ptr< Solver< SubstituteOperatorType > > substituteSolver,
                                   SubstituteOperatorType&                             substituteOperator )
   : holdsReference_( true )
   , substituteSolver_( substituteSolver )
   , substituteOperatorRef_( substituteOperator )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_UNUSED( A );
      if ( holdsReference_ )
      {
         for ( uint_t k = 0; k < x.getDimension(); k++ )
         {
            substituteSolver_->solve( substituteOperatorRef_.get(), x[k], b[k], level );
         }
      }
      else
      {
         for ( uint_t k = 0; k < x.getDimension(); k++ )
         {
            substituteSolver_->solve( *substituteOperator_, x[k], b[k], level );
         }
      }
   }

 private:
   const bool                                          holdsReference_;
   std::shared_ptr< Solver< SubstituteOperatorType > > substituteSolver_;
   std::shared_ptr< SubstituteOperatorType >           substituteOperator_;
   std::reference_wrapper< SubstituteOperatorType >    substituteOperatorRef_;
};

template < class OperatorType >
class ApplyPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   ApplyPreconditioner( std::shared_ptr< OperatorType > substituteOperator,
                        hyteg::DoFType                  flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   : holdsReference_( false )
   , substituteOperator_( substituteOperator )
   , substituteOperatorRef_( *substituteOperator )
   , flag_( flag )
   {}

   ApplyPreconditioner( OperatorType&  substituteOperator,
                        hyteg::DoFType flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   : holdsReference_( true )
   , substituteOperatorRef_( substituteOperator )
   , flag_( flag )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_UNUSED( A );
      if ( holdsReference_ )
      {
         substituteOperatorRef_.get().apply( b, x, level, flag_ );
      }
      else
      {
         substituteOperator_->apply( b, x, level, flag_ );
      }
   }

 private:
   const bool                             holdsReference_;
   std::shared_ptr< OperatorType >        substituteOperator_;
   std::reference_wrapper< OperatorType > substituteOperatorRef_;
   hyteg::DoFType                         flag_;
};

template < class OperatorType, class SubstituteOperatorTypeVelocity, class SubstituteOperatorTypePressure >
class TaylorHoodComponentwiseSubstitutePreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   TaylorHoodComponentwiseSubstitutePreconditioner(
       std::shared_ptr< Solver< SubstituteOperatorTypeVelocity > > substituteSolverVelocity,
       std::shared_ptr< SubstituteOperatorTypeVelocity >           substituteOperatorVelocity,
       std::shared_ptr< Solver< SubstituteOperatorTypePressure > > substituteSolverPressure   = nullptr,
       std::shared_ptr< SubstituteOperatorTypePressure >           substituteOperatorPressure = nullptr,
       bool                                                        projectMean                = true,
       std::shared_ptr< P2ProjectNormalOperator >                  projection                 = nullptr,
       DoFType                                                     projectionFlag             = FreeslipBoundary )
   : substituteSolverVelocity_( substituteSolverVelocity )
   , substituteOperatorVelocity_( substituteOperatorVelocity )
   , substituteSolverPressure_( substituteSolverPressure )
   , substituteOperatorPressure_( substituteOperatorPressure )
   , projectMean_( projectMean )
   , projectionFlag_( projectionFlag )
   , projection_( projection )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_UNUSED( A );
      if ( substituteSolverVelocity_ != nullptr && substituteOperatorVelocity_ != nullptr )
      {
         for ( uint_t k = 0; k < x.uvw().getDimension(); k++ )
         {
            substituteSolverVelocity_->solve( *substituteOperatorVelocity_, x.uvw()[k], b.uvw()[k], level );
         }
      }
      if ( substituteSolverPressure_ != nullptr && substituteOperatorPressure_ != nullptr )
      {
         substituteSolverPressure_->solve( *substituteOperatorPressure_, x.p(), b.p(), level );
      }

      if ( projectMean_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      if ( projection_ != nullptr )
      {
         projection_->project( x.uvw(), level, projectionFlag_ );
      }
   }

 private:
   std::shared_ptr< Solver< SubstituteOperatorTypeVelocity > > substituteSolverVelocity_;
   std::shared_ptr< SubstituteOperatorTypeVelocity >           substituteOperatorVelocity_;
   std::shared_ptr< Solver< SubstituteOperatorTypePressure > > substituteSolverPressure_;
   std::shared_ptr< SubstituteOperatorTypePressure >           substituteOperatorPressure_;

   bool    projectMean_;
   DoFType projectionFlag_;

   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

template < class OperatorType, class SubstituteOperatorTypeVelocity, class SubstituteOperatorTypePressure >
class TaylorHoodVectorSubstitutePreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   TaylorHoodVectorSubstitutePreconditioner(
       std::shared_ptr< Solver< SubstituteOperatorTypeVelocity > > substituteSolverVelocity,
       std::shared_ptr< SubstituteOperatorTypeVelocity >           substituteOperatorVelocity,
       std::shared_ptr< Solver< SubstituteOperatorTypePressure > > substituteSolverPressure   = nullptr,
       std::shared_ptr< SubstituteOperatorTypePressure >           substituteOperatorPressure = nullptr,
       bool                                                        projectMean                = true,
       std::shared_ptr< P2ProjectNormalOperator >                  projection                 = nullptr,
       DoFType                                                     projectionFlag             = FreeslipBoundary )
   : substituteSolverVelocity_( substituteSolverVelocity )
   , substituteOperatorVelocity_( substituteOperatorVelocity )
   , substituteSolverPressure_( substituteSolverPressure )
   , substituteOperatorPressure_( substituteOperatorPressure )
   , projectMean_( projectMean )
   , projectionFlag_( projectionFlag )
   , projection_( projection )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_UNUSED( A );
      if ( substituteSolverVelocity_ != nullptr && substituteOperatorVelocity_ != nullptr )
      {
         substituteSolverVelocity_->solve( *substituteOperatorVelocity_, x.uvw(), b.uvw(), level );
      }
      if ( substituteSolverPressure_ != nullptr && substituteOperatorPressure_ != nullptr )
      {
         substituteSolverPressure_->solve( *substituteOperatorPressure_, x.p(), b.p(), level );
      }

      if ( projectMean_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      if ( projection_ != nullptr )
      {
         projection_->project( x.uvw(), level, projectionFlag_ );
      }
   }

 private:
   std::shared_ptr< Solver< SubstituteOperatorTypeVelocity > > substituteSolverVelocity_;
   std::shared_ptr< SubstituteOperatorTypeVelocity >           substituteOperatorVelocity_;
   std::shared_ptr< Solver< SubstituteOperatorTypePressure > > substituteSolverPressure_;
   std::shared_ptr< SubstituteOperatorTypePressure >           substituteOperatorPressure_;

   bool    projectMean_;
   DoFType projectionFlag_;

   std::shared_ptr< P2ProjectNormalOperator > projection_;
};

} // namespace hyteg