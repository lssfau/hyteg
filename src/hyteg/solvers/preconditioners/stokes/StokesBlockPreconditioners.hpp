/*
 * Copyright (c) 2024 Andreas Burkhart.
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

#include "core/math/Random.h"

#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/ZeroOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hyteg {

// -------------------------------------------------------------
// -------------------------- Preface --------------------------
// -------------------------------------------------------------

// The block preconditioners in this header expect an operator of the form
//    {{A, B^T},{B, 0}}
//
// The operator in question is expected to operate on functions that provide
// access to their respective parts via *.uvw() and *.p() with typedefs
// VelocityFunction_T and PressureFunction_T for the respective types.
//
// Additionally the operator is expected to define the typedefs
//    AOperatorType             for the ABlock                 Operator,
//    BOperatorType             for the BBlock                 Operator,
//    BTOperatorType            for the BTBlock                Operator,
//    SchurOperatorType         for the Schur complement       Operator and
//    StabilisationOperatorType for the optional stabilisation Operator
//
// and provide access to its inner operators via the functions
//    const AOperatorType&              getA()     const
//    const BOperatorType&              getB()     const
//    const BTOperatorType&             getBT()    const
//    const StabilisationOperatorType&  getStab()  const
//
// If you do not want to use solvers für the approximations but apply operators instead
// (e.g. invLumpedMass, PSPG) you can use the the
//    ApplyPreconditioner<OperatorType>
// class found in
//    src/hyteg/solvers/SubstitutePreconditioner.hpp
//
// If you do not want to use stabilisation, you can substitute the operator
// with the
//    ZeroOperator< typename SrcType, typename DstType>
// class found in src/hyteg/operators/.

// --------------------------------------------------------------
// ---------------------------- Note ----------------------------
// --------------------------------------------------------------

// The loops applying multiple velocity/pressure updated could also be achieved via loops of the form
//
// A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
// residual_.uvw().assign( { 1.0, -1.0 }, { b.uvw(), residual_.uvw() }, level, flag_ );
//
// for ( uint_t i = 0; i < 2; i++ )
// {
//    ABlockApproximationSolver_->solve( A.getA(), x.uvw(), residual_.uvw(), level );
// }
//
// but here we want to specifically enforce e.g. the usage of A.getA() for the residual calculation
// (some solver wrappers just ignore the given operator and use something else).

// -------------------------------------------------------------
// ---------------------- Preconditioners ----------------------
// -------------------------------------------------------------

template < typename OperatorType,
           typename AOperatorType,
           typename SchurOperatorType,
           class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class InexactUzawaPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   InexactUzawaPreconditioner( const std::shared_ptr< PrimitiveStorage >&            storage,
                               const uint_t                                          minLevel,
                               const uint_t                                          maxLevel,
                               const SchurOperatorType&                              schurOp,
                               const std::shared_ptr< Solver< AOperatorType > >&     ABlockApproximationSolver,
                               const std::shared_ptr< Solver< SchurOperatorType > >& SchurComplementApproximationSolver,
                               real_t                                                relaxParamA,
                               real_t                                                relaxParamSchur,
                               uint_t                                                VelocityIterations = 1,
                               std::shared_ptr< VelocityProjectionOperatorType >     projection         = nullptr,
                               hyteg::DoFType flag            = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
                               bool           projectPressure = true )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( schurOp )
   , ABlockApproximationSolver_( ABlockApproximationSolver )
   , SchurComplementApproximationSolver_( SchurComplementApproximationSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , projectPressure_( projectPressure )
   , residual_( FunctionType( "InexactUzawaPreconditioner residual", storage, minLevel, maxLevel ) )
   , tmp_( FunctionType( "InexactUzawaPreconditioner tmp", storage, minLevel, maxLevel ) )
   , projection_( projection )
   {}

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      residual_.copyBoundaryConditionFromFunction( x );
      tmp_.copyBoundaryConditionFromFunction( x );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      if ( projection_ != nullptr )
      {
         projection_->project( x.uvw(), level, FreeslipBoundary );
      }

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual
         A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
         A.getA().apply( x.uvw(), residual_.uvw(), level, flag_, Add );
         residual_.uvw().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.uvw(), residual_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( residual_.uvw(), level, FreeslipBoundary );
         }

         // apply A_hat
         tmp_.uvw().interpolate( real_c( 0.0 ), level, flag_ );
         ABlockApproximationSolver_->solve( A.getA(), tmp_.uvw(), residual_.uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1.0 ), relaxParamA_ }, { x.uvw(), tmp_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( x.uvw(), level, FreeslipBoundary );
         }
      }

      // --------------------------Pressure--------------------------

      // calculate residual
      A.getB().apply( x.uvw(), residual_.p(), level, flag_, Replace );
      A.getStab().apply( x.p(), residual_.p(), level, flag_, Add );
      residual_.p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residual_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( residual_.p(), level );
      }

      // apply S_hat
      tmp_.p().interpolate( real_c( 0.0 ), level, flag_ );
      SchurComplementApproximationSolver_->solve( schurOp_, tmp_.p(), residual_.p(), level );

      // update x with the scaled result
      x.p().assign( { real_c( 1.0 ), -relaxParamSchur_ }, { x.p(), tmp_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }
   }

 protected:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;

   SchurOperatorType schurOp_;

   std::shared_ptr< Solver< AOperatorType > >     ABlockApproximationSolver_;
   std::shared_ptr< Solver< SchurOperatorType > > SchurComplementApproximationSolver_;

   DoFType flag_;
   bool    hasGlobalCells_;

   real_t       relaxParamA_;
   real_t       relaxParamSchur_;
   uint_t       VelocityIterations_;
   bool         projectPressure_;
   FunctionType residual_;
   FunctionType tmp_;

   std::shared_ptr< VelocityProjectionOperatorType > projection_;
};

template < typename OperatorType,
           typename AOperatorType,
           typename SchurOperatorType,
           class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class AdjointInexactUzawaPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   AdjointInexactUzawaPreconditioner( const std::shared_ptr< PrimitiveStorage >&            storage,
                                      const uint_t                                          minLevel,
                                      const uint_t                                          maxLevel,
                                      const SchurOperatorType&                              schurOp,
                                      const std::shared_ptr< Solver< AOperatorType > >&     ABlockApproximationSolver,
                                      const std::shared_ptr< Solver< SchurOperatorType > >& SchurComplementApproximationSolver,
                                      real_t                                                relaxParamA,
                                      real_t                                                relaxParamSchur,
                                      uint_t                                                VelocityIterations = 1,
                                      std::shared_ptr< VelocityProjectionOperatorType >     projection         = nullptr,
                                      hyteg::DoFType flag = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
                                      bool           projectPressure = true )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( schurOp )
   , ABlockApproximationSolver_( ABlockApproximationSolver )
   , SchurComplementApproximationSolver_( SchurComplementApproximationSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , projectPressure_( projectPressure )
   , residual_( FunctionType( "AdjointInexactUzawaPreconditioner residual", storage, minLevel, maxLevel ) )
   , tmp_( FunctionType( "AdjointInexactUzawaPreconditioner tmp", storage, minLevel, maxLevel ) )
   , projection_( projection )
   {}

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      residual_.copyBoundaryConditionFromFunction( x );
      tmp_.copyBoundaryConditionFromFunction( x );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      if ( projection_ != nullptr )
      {
         projection_->project( x.uvw(), level, FreeslipBoundary );
      }

      // --------------------------Pressure--------------------------

      // calculate residual
      A.getB().apply( x.uvw(), residual_.p(), level, flag_, Replace );
      A.getStab().apply( x.p(), residual_.p(), level, flag_, Add );
      residual_.p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residual_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( residual_.p(), level );
      }

      // apply S_hat
      tmp_.p().interpolate( real_c( 0.0 ), level, flag_ );
      SchurComplementApproximationSolver_->solve( schurOp_, tmp_.p(), residual_.p(), level );

      // update x with the scaled result
      x.p().assign( { real_c( 1.0 ), -relaxParamSchur_ }, { x.p(), tmp_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      // --------------------------Velocity--------------------------

      if ( projection_ != nullptr )
      {
         projection_->project( x.uvw(), level, FreeslipBoundary );
      }

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual
         A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
         A.getA().apply( x.uvw(), residual_.uvw(), level, flag_, Add );
         residual_.uvw().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.uvw(), residual_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( residual_.uvw(), level, FreeslipBoundary );
         }

         // apply A_hat
         tmp_.uvw().interpolate( real_c( 0.0 ), level, flag_ );
         ABlockApproximationSolver_->solve( A.getA(), tmp_.uvw(), residual_.uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1.0 ), relaxParamA_ }, { x.uvw(), tmp_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( x.uvw(), level, FreeslipBoundary );
         }
      }
   }

 protected:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;

   SchurOperatorType schurOp_;

   std::shared_ptr< Solver< AOperatorType > >     ABlockApproximationSolver_;
   std::shared_ptr< Solver< SchurOperatorType > > SchurComplementApproximationSolver_;

   DoFType flag_;
   bool    hasGlobalCells_;

   real_t       relaxParamA_;
   real_t       relaxParamSchur_;
   uint_t       VelocityIterations_;
   bool         projectPressure_;
   FunctionType residual_;
   FunctionType tmp_;

   std::shared_ptr< VelocityProjectionOperatorType > projection_;
};

template < typename OperatorType,
           typename AOperatorType,
           typename SchurOperatorType,
           class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class BlockApproximateFactorisationPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   BlockApproximateFactorisationPreconditioner(
       const std::shared_ptr< PrimitiveStorage >&            storage,
       const uint_t                                          minLevel,
       const uint_t                                          maxLevel,
       const SchurOperatorType&                              schurOp,
       const std::shared_ptr< Solver< AOperatorType > >&     ABlockApproximationSolver,
       const std::shared_ptr< Solver< SchurOperatorType > >& SchurComplementApproximationSolver,
       real_t                                                relaxParamA,
       real_t                                                relaxParamSchur,
       uint_t                                                VelocityIterations = 1,
       std::shared_ptr< VelocityProjectionOperatorType >     projection         = nullptr,
       hyteg::DoFType flag            = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
       bool           projectPressure = true )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( schurOp )
   , ABlockApproximationSolver_( ABlockApproximationSolver )
   , SchurComplementApproximationSolver_( SchurComplementApproximationSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , projectPressure_( projectPressure )
   , residual_( FunctionType( "BlockApproximateFactorisationPreconditioner residual", storage, minLevel, maxLevel ) )
   , tmp_( FunctionType( "BlockApproximateFactorisationPreconditioner tmp", storage, minLevel, maxLevel ) )
   , projection_( projection )
   {}

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      residual_.copyBoundaryConditionFromFunction( x );
      tmp_.copyBoundaryConditionFromFunction( x );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      if ( projection_ != nullptr )
      {
         projection_->project( x.uvw(), level, FreeslipBoundary );
      }

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual
         A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
         A.getA().apply( x.uvw(), residual_.uvw(), level, flag_, Add );
         residual_.uvw().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.uvw(), residual_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( residual_.uvw(), level, FreeslipBoundary );
         }

         // apply A_hat
         tmp_.uvw().interpolate( real_c( 0.0 ), level, flag_ );
         ABlockApproximationSolver_->solve( A.getA(), tmp_.uvw(), residual_.uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1.0 ), relaxParamA_ }, { x.uvw(), tmp_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( x.uvw(), level, FreeslipBoundary );
         }
      }

      // --------------------------Pressure--------------------------

      // calculate residual
      A.getB().apply( x.uvw(), residual_.p(), level, flag_, Replace );
      A.getStab().apply( x.p(), residual_.p(), level, flag_, Add );
      residual_.p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residual_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( residual_.p(), level );
      }

      // apply S_hat
      tmp_.p().interpolate( real_c( 0.0 ), level, flag_ );
      SchurComplementApproximationSolver_->solve( schurOp_, tmp_.p(), residual_.p(), level );
      tmp_.p().assign( { -relaxParamSchur_ }, { tmp_.p() }, level, flag_ );

      // update x with the scaled result
      x.p().assign( { real_c( 1.0 ), real_c( 1.0 ) }, { x.p(), tmp_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      // --------------------------Velocity--------------------------

      {
         if ( projectPressure_ )
         {
            vertexdof::projectMean( tmp_.p(), level );
         }

         A.getBT().apply( tmp_.p(), residual_.uvw(), level, flag_, Replace );

         if ( projection_ != nullptr )
         {
            projection_->project( residual_.uvw(), level, FreeslipBoundary );
         }

         // apply A_hat
         tmp_.uvw().interpolate( real_c( 0.0 ), level, flag_ );
         ABlockApproximationSolver_->solve( A.getA(), tmp_.uvw(), residual_.uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1.0 ), -relaxParamA_ }, { x.uvw(), tmp_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( x.uvw(), level, FreeslipBoundary );
         }
      }

      for ( uint_t i = 0; i < VelocityIterations_ - 1; i++ )
      {
         // calculate residual
         A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
         A.getA().apply( x.uvw(), residual_.uvw(), level, flag_, Add );
         residual_.uvw().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.uvw(), residual_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( residual_.uvw(), level, FreeslipBoundary );
         }

         // apply A_hat
         tmp_.uvw().interpolate( real_c( 0.0 ), level, flag_ );
         ABlockApproximationSolver_->solve( A.getA(), tmp_.uvw(), residual_.uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1.0 ), relaxParamA_ }, { x.uvw(), tmp_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( x.uvw(), level, FreeslipBoundary );
         }
      }
   }

 protected:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;

   SchurOperatorType schurOp_;

   std::shared_ptr< Solver< AOperatorType > >     ABlockApproximationSolver_;
   std::shared_ptr< Solver< SchurOperatorType > > SchurComplementApproximationSolver_;

   DoFType flag_;
   bool    hasGlobalCells_;

   real_t       relaxParamA_;
   real_t       relaxParamSchur_;
   uint_t       VelocityIterations_;
   bool         projectPressure_;
   FunctionType residual_;
   FunctionType tmp_;

   std::shared_ptr< VelocityProjectionOperatorType > projection_;
};

template < typename OperatorType,
           typename AOperatorType,
           typename SchurOperatorType,
           class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class SymmetricUzawaPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   SymmetricUzawaPreconditioner( const std::shared_ptr< PrimitiveStorage >&            storage,
                                 const uint_t                                          minLevel,
                                 const uint_t                                          maxLevel,
                                 const SchurOperatorType&                              schurOp,
                                 const std::shared_ptr< Solver< AOperatorType > >&     ABlockApproximationSolver,
                                 const std::shared_ptr< Solver< SchurOperatorType > >& SchurComplementApproximationSolver,
                                 real_t                                                relaxParamA,
                                 real_t                                                relaxParamSchur,
                                 uint_t                                                VelocityIterations = 1,
                                 std::shared_ptr< VelocityProjectionOperatorType >     projection         = nullptr,
                                 hyteg::DoFType flag            = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
                                 bool           projectPressure = true )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , schurOp_( schurOp )
   , ABlockApproximationSolver_( ABlockApproximationSolver )
   , SchurComplementApproximationSolver_( SchurComplementApproximationSolver )
   , flag_( flag )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParamA_( relaxParamA )
   , relaxParamSchur_( relaxParamSchur )
   , VelocityIterations_( VelocityIterations )
   , projectPressure_( projectPressure )
   , residual_( FunctionType( "SymmetricUzawaPreconditioner residual", storage, minLevel, maxLevel ) )
   , tmp_( FunctionType( "SymmetricUzawaPreconditioner tmp", storage, minLevel, maxLevel ) )
   , projection_( projection )
   {}

   virtual void solve( const OperatorType&                   A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const uint_t                          level ) override
   {
      residual_.copyBoundaryConditionFromFunction( x );
      tmp_.copyBoundaryConditionFromFunction( x );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      if ( projection_ != nullptr )
      {
         projection_->project( x.uvw(), level, FreeslipBoundary );
      }

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual
         A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
         A.getA().apply( x.uvw(), residual_.uvw(), level, flag_, Add );
         residual_.uvw().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.uvw(), residual_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( residual_.uvw(), level, FreeslipBoundary );
         }

         // apply A_hat
         tmp_.uvw().interpolate( real_c( 0.0 ), level, flag_ );
         ABlockApproximationSolver_->solve( A.getA(), tmp_.uvw(), residual_.uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1.0 ), relaxParamA_ }, { x.uvw(), tmp_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( x.uvw(), level, FreeslipBoundary );
         }
      }

      // --------------------------Pressure--------------------------

      // calculate residual
      A.getB().apply( x.uvw(), residual_.p(), level, flag_, Replace );
      A.getStab().apply( x.p(), residual_.p(), level, flag_, Add );
      residual_.p().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.p(), residual_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( residual_.p(), level );
      }

      // apply S_hat
      tmp_.p().interpolate( real_c( 0.0 ), level, flag_ );
      SchurComplementApproximationSolver_->solve( schurOp_, tmp_.p(), residual_.p(), level );

      // update x with the scaled result
      x.p().assign( { real_c( 1.0 ), -relaxParamSchur_ }, { x.p(), tmp_.p() }, level, flag_ );

      if ( projectPressure_ )
      {
         vertexdof::projectMean( x.p(), level );
      }

      // --------------------------Velocity--------------------------

      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         // calculate residual
         A.getBT().apply( x.p(), residual_.uvw(), level, flag_, Replace );
         A.getA().apply( x.uvw(), residual_.uvw(), level, flag_, Add );
         residual_.uvw().assign( { real_c( 1.0 ), real_c( -1.0 ) }, { b.uvw(), residual_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( residual_.uvw(), level, FreeslipBoundary );
         }

         // apply A_hat
         tmp_.uvw().interpolate( real_c( 0.0 ), level, flag_ );
         ABlockApproximationSolver_->solve( A.getA(), tmp_.uvw(), residual_.uvw(), level );

         // update x with the scaled result
         x.uvw().assign( { real_c( 1.0 ), relaxParamA_ }, { x.uvw(), tmp_.uvw() }, level, flag_ );

         if ( projection_ != nullptr )
         {
            projection_->project( x.uvw(), level, FreeslipBoundary );
         }
      }
   }

 protected:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;

   SchurOperatorType schurOp_;

   std::shared_ptr< Solver< AOperatorType > >     ABlockApproximationSolver_;
   std::shared_ptr< Solver< SchurOperatorType > > SchurComplementApproximationSolver_;

   DoFType flag_;
   bool    hasGlobalCells_;

   real_t       relaxParamA_;
   real_t       relaxParamSchur_;
   uint_t       VelocityIterations_;
   bool         projectPressure_;
   FunctionType residual_;
   FunctionType tmp_;

   std::shared_ptr< VelocityProjectionOperatorType > projection_;
};

// -------------------------------------------------------------
// -------------------------- Utility --------------------------
// -------------------------------------------------------------

template < typename OperatorType,
           typename AOperatorType,
           typename SchurOperatorType,
           class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class UzawaOmegaEstimationOperator
: public Operator< typename OperatorType::srcType::PressureFunction_T, typename OperatorType::dstType::PressureFunction_T >
{
 public:
   UzawaOmegaEstimationOperator( const std::shared_ptr< PrimitiveStorage >&            storage,
                                 uint_t                                                minLevel,
                                 uint_t                                                maxLevel,
                                 const OperatorType&                                   StokesOp,
                                 const SchurOperatorType&                              schurOp,
                                 const std::shared_ptr< Solver< AOperatorType > >&     ABlockSolver,
                                 const std::shared_ptr< Solver< SchurOperatorType > >& SchurSolver,
                                 const uint_t                                          VelocityIterations = 1,
                                 BoundaryCondition                                     bc                 = BoundaryCondition(),
                                 std::shared_ptr< VelocityProjectionOperatorType >     projection         = nullptr )
   : UzawaOmegaEstimationOperator( storage,
                                   minLevel,
                                   maxLevel,
                                   StokesOp,
                                   schurOp,
                                   ABlockSolver,
                                   SchurSolver,
                                   VelocityIterations,
                                   bc,
                                   bc,
                                   bc,
                                   projection )
   {}

   UzawaOmegaEstimationOperator( const std::shared_ptr< PrimitiveStorage >&            storage,
                                 uint_t                                                minLevel,
                                 uint_t                                                maxLevel,
                                 const OperatorType&                                   StokesOp,
                                 const SchurOperatorType&                              schurOp,
                                 const std::shared_ptr< Solver< AOperatorType > >&     ABlockSolver,
                                 const std::shared_ptr< Solver< SchurOperatorType > >& SchurSolver,
                                 const uint_t                                          VelocityIterations = 1,
                                 BoundaryCondition                                     bcX                = BoundaryCondition(),
                                 BoundaryCondition                                     bcY                = BoundaryCondition(),
                                 BoundaryCondition                                     bcZ                = BoundaryCondition(),
                                 std::shared_ptr< VelocityProjectionOperatorType >     projection         = nullptr )
   : Operator< typename OperatorType::srcType::PressureFunction_T, typename OperatorType::dstType::PressureFunction_T >(
         storage,
         minLevel,
         maxLevel )
   , StokesOp_( StokesOp )
   , schurOp_( schurOp )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , ABlockSolver_( ABlockSolver )
   , SchurSolver_( SchurSolver )
   , VelocityIterations_( VelocityIterations )
   , tmp_rhs_( "tmp_rhs_", storage, minLevel, maxLevel )
   , tmp_solution_( "tmp_solution_", storage, minLevel, maxLevel )
   , tmp_schur_( "tmp_schur", storage, minLevel, maxLevel )
   , projection_( projection )
   {
      tmp_rhs_[0].setBoundaryCondition( bcX );
      tmp_rhs_[1].setBoundaryCondition( bcY );
      if ( storage->hasGlobalCells() )
      {
         tmp_rhs_[2].setBoundaryCondition( bcZ );
      }

      tmp_solution_[0].setBoundaryCondition( bcX );
      tmp_solution_[1].setBoundaryCondition( bcY );
      if ( storage->hasGlobalCells() )
      {
         tmp_solution_[2].setBoundaryCondition( bcZ );
      }
   }

   void apply( const typename OperatorType::srcType::PressureFunction_T& src,
               const typename OperatorType::dstType::PressureFunction_T& dst,
               const uint_t                                              level,
               const DoFType                                             flag,
               UpdateType                                                updateType = Replace ) const
   {
      WALBERLA_UNUSED( updateType );
      vertexdof::projectMean( src, level );

      StokesOp_.getBT().apply( src, tmp_rhs_, level, flag, Replace );

      if ( projection_ != nullptr )
      {
         projection_->project( tmp_rhs_, level, FreeslipBoundary );
      }

      tmp_solution_.interpolate( real_c( 0.0 ), level, All );
      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         ABlockSolver_->solve( StokesOp_.getA(), tmp_solution_, tmp_rhs_, level );

         if ( projection_ != nullptr )
         {
            projection_->project( tmp_solution_, level, FreeslipBoundary );
         }
      }

      StokesOp_.getStab().apply( src, tmp_schur_, level, flag, Replace );
      tmp_schur_.assign( { real_c( -1.0 ) }, { tmp_schur_ }, level, flag );

      StokesOp_.getB().apply( tmp_solution_, tmp_schur_, level, flag, Add );

      vertexdof::projectMean( tmp_schur_, level );

      dst.interpolate( real_c( 0.0 ), level, flag );
      SchurSolver_->solve( schurOp_, dst, tmp_schur_, level );

      vertexdof::projectMean( dst, level );
   }

   OperatorType      StokesOp_;
   SchurOperatorType schurOp_;

   bool                                           hasGlobalCells_;
   std::shared_ptr< Solver< AOperatorType > >     ABlockSolver_;
   std::shared_ptr< Solver< SchurOperatorType > > SchurSolver_;

   uint_t VelocityIterations_;

   typename OperatorType::srcType::VelocityFunction_T tmp_rhs_;
   typename OperatorType::srcType::VelocityFunction_T tmp_solution_;
   typename OperatorType::srcType::PressureFunction_T tmp_schur_;

   std::shared_ptr< VelocityProjectionOperatorType > projection_;
};

// This function tries to calculate a scaling for the Schur complement approximation S_hat, s.t. S_hat >= S is fulfilled
template < typename OperatorType,
           typename AOperatorType,
           typename SchurOperatorType,
           class VelocityProjectionOperatorType = P2ProjectNormalOperator >
real_t estimateUzawaOmega( const std::shared_ptr< PrimitiveStorage >&            storage,
                           const uint_t&                                         minLevel,
                           const uint_t&                                         maxLevel,
                           const OperatorType&                                   StokesOp,
                           const SchurOperatorType&                              schurOp,
                           const std::shared_ptr< Solver< AOperatorType > >&     ABlockSolver,
                           const std::shared_ptr< Solver< SchurOperatorType > >& SchurSolver,
                           const uint_t&                                         numPowerIterations,
                           const uint_t                                          VelocityIterations = 1,
                           BoundaryCondition                                     bcX                = BoundaryCondition(),
                           BoundaryCondition                                     bcY                = BoundaryCondition(),
                           BoundaryCondition                                     bcZ                = BoundaryCondition(),
                           uint_fast32_t                                         randomSeed         = 42,
                           std::shared_ptr< VelocityProjectionOperatorType >     projection         = nullptr )
{
   UzawaOmegaEstimationOperator< OperatorType, AOperatorType, SchurOperatorType > estimator(
       storage, minLevel, maxLevel, StokesOp, schurOp, ABlockSolver, SchurSolver, VelocityIterations, bcX, bcY, bcZ, projection );
   typename OperatorType::srcType::PressureFunction_T iterationVector( "iterationVector", storage, minLevel, maxLevel );
   typename OperatorType::srcType::PressureFunction_T auxVector( "auxVector", storage, minLevel, maxLevel );
   walberla::math::seedRandomGenerator( randomSeed );
   auto randFunction = []( const Point3D& ) { return walberla::math::realRandom(); };
   iterationVector.interpolate( randFunction, maxLevel, All );
   const real_t estimatedRelaxationParameter =
       estimateSpectralRadiusWithPowerIteration( estimator, iterationVector, auxVector, numPowerIterations, storage, maxLevel );
   return real_c( 1.0 ) / estimatedRelaxationParameter;
}

template < typename OperatorType, typename AOperatorType, class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class UzawaSigmaEstimationOperator
: public Operator< typename OperatorType::srcType::VelocityFunction_T, typename OperatorType::dstType::VelocityFunction_T >
{
 public:
   UzawaSigmaEstimationOperator( const std::shared_ptr< PrimitiveStorage >&        storage,
                                 uint_t                                            minLevel,
                                 uint_t                                            maxLevel,
                                 const OperatorType&                               StokesOp,
                                 const std::shared_ptr< Solver< AOperatorType > >& ABlockSolver,
                                 const uint_t                                      VelocityIterations = 1,
                                 BoundaryCondition                                 bc                 = BoundaryCondition(),
                                 std::shared_ptr< VelocityProjectionOperatorType > projection         = nullptr )
   : UzawaSigmaEstimationOperator( storage,
                                   minLevel,
                                   maxLevel,
                                   StokesOp,
                                   ABlockSolver,
                                   VelocityIterations,
                                   bc,
                                   bc,
                                   bc,
                                   projection )
   {}

   UzawaSigmaEstimationOperator( const std::shared_ptr< PrimitiveStorage >&        storage,
                                 uint_t                                            minLevel,
                                 uint_t                                            maxLevel,
                                 const OperatorType&                               StokesOp,
                                 const std::shared_ptr< Solver< AOperatorType > >& ABlockSolver,
                                 const uint_t                                      VelocityIterations = 1,
                                 BoundaryCondition                                 bcX                = BoundaryCondition(),
                                 BoundaryCondition                                 bcY                = BoundaryCondition(),
                                 BoundaryCondition                                 bcZ                = BoundaryCondition(),
                                 std::shared_ptr< VelocityProjectionOperatorType > projection         = nullptr )
   : Operator< typename OperatorType::srcType::VelocityFunction_T, typename OperatorType::dstType::VelocityFunction_T >(
         storage,
         minLevel,
         maxLevel )
   , StokesOp_( StokesOp )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , ABlockSolver_( ABlockSolver )
   , VelocityIterations_( VelocityIterations )
   , tmp_rhs_( "tmp_rhs_", storage, minLevel, maxLevel )
   , projection_( projection )
   {
      tmp_rhs_[0].setBoundaryCondition( bcX );
      tmp_rhs_[1].setBoundaryCondition( bcY );
      if ( storage->hasGlobalCells() )
      {
         tmp_rhs_[2].setBoundaryCondition( bcZ );
      }
   }

   void apply( const typename OperatorType::srcType::VelocityFunction_T& src,
               const typename OperatorType::dstType::VelocityFunction_T& dst,
               const uint_t                                              level,
               const DoFType                                             flag,
               UpdateType                                                updateType = Replace ) const
   {
      WALBERLA_UNUSED( updateType );

      if ( projection_ != nullptr )
      {
         projection_->project( src, level, FreeslipBoundary );
      }

      StokesOp_.getA().apply( src, tmp_rhs_, level, flag, Replace );

      if ( projection_ != nullptr )
      {
         projection_->project( tmp_rhs_, level, FreeslipBoundary );
      }

      dst.interpolate( real_c( 0.0 ), level, All );
      for ( uint_t i = 0; i < VelocityIterations_; i++ )
      {
         ABlockSolver_->solve( StokesOp_.getA(), dst, tmp_rhs_, level );

         if ( projection_ != nullptr )
         {
            projection_->project( dst, level, FreeslipBoundary );
         }
      }
   }

   OperatorType StokesOp_;

   bool                                       hasGlobalCells_;
   std::shared_ptr< Solver< AOperatorType > > ABlockSolver_;

   uint_t VelocityIterations_;

   typename OperatorType::srcType::VelocityFunction_T tmp_rhs_;

   std::shared_ptr< VelocityProjectionOperatorType > projection_;
};

// This function tries to calculate a scaling for the A Block approximation A_hat, s.t. A_hat >= A is fulfilled
template < typename OperatorType, typename AOperatorType, class VelocityProjectionOperatorType = P2ProjectNormalOperator >
real_t estimateUzawaSigma( const std::shared_ptr< PrimitiveStorage >&        storage,
                           const uint_t&                                     minLevel,
                           const uint_t&                                     maxLevel,
                           const OperatorType&                               StokesOp,
                           const std::shared_ptr< Solver< AOperatorType > >& ABlockSolver,
                           const uint_t&                                     numPowerIterations,
                           const uint_t                                      VelocityIterations = 1,
                           BoundaryCondition                                 bcX                = BoundaryCondition(),
                           BoundaryCondition                                 bcY                = BoundaryCondition(),
                           BoundaryCondition                                 bcZ                = BoundaryCondition(),
                           uint_fast32_t                                     randomSeed         = 42,
                           std::shared_ptr< VelocityProjectionOperatorType > projection         = nullptr )
{
   UzawaSigmaEstimationOperator< OperatorType, AOperatorType > estimator(
       storage, minLevel, maxLevel, StokesOp, ABlockSolver, VelocityIterations, bcX, bcY, bcZ, projection );
   typename OperatorType::srcType::VelocityFunction_T iterationVector( "iterationVector", storage, minLevel, maxLevel );
   typename OperatorType::srcType::VelocityFunction_T auxVector( "auxVector", storage, minLevel, maxLevel );
   walberla::math::seedRandomGenerator( randomSeed );
   auto randFunction = []( const Point3D& ) { return walberla::math::realRandom(); };
   iterationVector.interpolate( randFunction, maxLevel, All );
   const real_t estimatedRelaxationParameter =
       estimateSpectralRadiusWithPowerIteration( estimator, iterationVector, auxVector, numPowerIterations, storage, maxLevel );
   return real_c( 1.0 ) / estimatedRelaxationParameter;
}

template < class OperatorType, class VelocityOperatorType >
class VelocityBlockWrapper : public Solver< OperatorType >
{
 public:
   VelocityBlockWrapper( const std::shared_ptr< PrimitiveStorage >&               storage,
                         const std::shared_ptr< Solver< VelocityOperatorType > >& velocitySolver,
                         const std::shared_ptr< VelocityOperatorType >&           velocityOperator )
   : hasGlobalCells_( storage->hasGlobalCells() )
   , velocitySolver_( velocitySolver )
   , velocityOperator_( velocityOperator )
   {}

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const uint_t                          level ) override
   {
      velocitySolver_->solve( *velocityOperator_, x.uvw(), b.uvw(), level );
   }

 private:
   bool                                              hasGlobalCells_;
   std::shared_ptr< Solver< VelocityOperatorType > > velocitySolver_;
   std::shared_ptr< VelocityOperatorType >           velocityOperator_;
};

template < typename OperatorType, typename AOperatorType, class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class SchurComplementOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   SchurComplementOperator( const std::shared_ptr< PrimitiveStorage >&        storage,
                            size_t                                            minLevel,
                            size_t                                            maxLevel,
                            const OperatorType&                               StokesOp,
                            const std::shared_ptr< Solver< AOperatorType > >& ABlockSolver,
                            BoundaryCondition                                 VelocityBC,
                            std::shared_ptr< VelocityProjectionOperatorType > projection = nullptr )
   : SchurComplementOperator( storage,
                              minLevel,
                              maxLevel,
                              StokesOp,
                              ABlockSolver,
                              VelocityBC,
                              VelocityBC,
                              VelocityBC,
                              projection )
   {}

   SchurComplementOperator( const std::shared_ptr< PrimitiveStorage >&        storage,
                            size_t                                            minLevel,
                            size_t                                            maxLevel,
                            const OperatorType&                               StokesOp,
                            const std::shared_ptr< Solver< AOperatorType > >& ABlockSolver,
                            BoundaryCondition                                 VelocityBCx,
                            BoundaryCondition                                 VelocityBCy,
                            BoundaryCondition                                 VelocityBCz,
                            std::shared_ptr< VelocityProjectionOperatorType > projection = nullptr )
   : Operator( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , temporary_( "schur temporary", storage, minLevel, maxLevel )
   , temporarySolution_( "schur temporary solution", storage, minLevel, maxLevel )
   , StokesOp_( StokesOp )
   , ABlockSolver_( ABlockSolver )
   , projection_( projection )
   {
      temporary_[0].setBoundaryCondition( VelocityBCx );
      temporary_[1].setBoundaryCondition( VelocityBCy );
      if ( hasGlobalCells_ )
      {
         temporary_[2].setBoundaryCondition( VelocityBCz );
      }

      temporarySolution_[0].setBoundaryCondition( VelocityBCx );
      temporarySolution_[1].setBoundaryCondition( VelocityBCy );
      if ( hasGlobalCells_ )
      {
         temporarySolution_[2].setBoundaryCondition( VelocityBCz );
      }
   }

   void apply( const typename OperatorType::srcType::PressureFunction_T& src,
               const typename OperatorType::srcType::PressureFunction_T& dst,
               const uint_t                                              level,
               const DoFType                                             flag,
               const UpdateType                                          updateType = Replace ) const
   {
      vertexdof::projectMean( src, level );

      StokesOp_.getBT().apply( src, temporary_, level, flag, Replace );

      if ( projection_ != nullptr )
      {
         projection_->project( temporary_, level, FreeslipBoundary );
      }

      temporarySolution_.interpolate( real_c( 0.0 ), level, flag );
      ABlockSolver_->solve( StokesOp_.getA(), temporarySolution_, temporary_, level );

      if ( projection_ != nullptr )
      {
         projection_->project( temporarySolution_, level, FreeslipBoundary );
      }

      if ( updateType == Replace )
      {
         StokesOp_.getB().apply( temporarySolution_, dst, level, flag, Replace );
         StokesOp_.getStab().apply( src, dst, level, flag, Add );
      }
      else
      {
         StokesOp_.getB().apply( temporarySolution_, dst, level, flag, Add );
         StokesOp_.getStab().apply( src, dst, level, flag, Add );
      }

      vertexdof::projectMean( dst, level );
   }

 private:
   bool                                                       hasGlobalCells_;
   mutable typename OperatorType::srcType::VelocityFunction_T temporary_;
   mutable typename OperatorType::srcType::VelocityFunction_T temporarySolution_;
   OperatorType                                               StokesOp_;
   std::shared_ptr< Solver< AOperatorType > >                 ABlockSolver_;

   std::shared_ptr< VelocityProjectionOperatorType > projection_;
};

template < typename OperatorType, typename SchurOperatorType, class VelocityProjectionOperatorType = P2ProjectNormalOperator >
class InverseSchurComplementPreconditioner : public Solver< typename OperatorType::SchurOperatorType >
{
 public:
   typedef typename SchurOperatorType::srcType FunctionType;

   InverseSchurComplementPreconditioner(
       const std::shared_ptr< PrimitiveStorage >&                                                                  storage,
       size_t                                                                                                      minLevel,
       size_t                                                                                                      maxLevel,
       const SchurComplementOperator< OperatorType, VelocityProjectionOperatorType >&                              SchurOperator,
       const std::shared_ptr< Solver< SchurComplementOperator< OperatorType, VelocityProjectionOperatorType > > >& SchurSolver )
   : SchurOperator_( SchurOperator )
   , SchurSolver_( SchurSolver )
   {}

   void solve( const SchurOperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_UNUSED( A );

      x.interpolate( real_c( 0.0 ), level, All );
      SchurSolver_->solve( SchurOperator_, x, b, level );
   }

 private:
   SchurComplementOperator< OperatorType, VelocityProjectionOperatorType >                              SchurOperator_;
   std::shared_ptr< Solver< SchurComplementOperator< OperatorType, VelocityProjectionOperatorType > > > SchurSolver_;
};

} // namespace hyteg
