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
#include "core/Environment.h"
#include "core/config/Config.h"

#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/misc/SFINAE.hpp"
#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FunctionMultiplicationPreconditioner.hpp"
#include "hyteg/solvers/Solver.hpp"

#include "../../OperatorTools/BFBT.hpp"
#include "../../OperatorTools/SchurOperator.hpp"
#include "SchurIdentityPreconditioner.hpp"
#include "SchurSolver.hpp"

namespace MantleConvection {

template < class OperatorType,
           class AApproximationOperatorTypeRight = typename OperatorType::AOperatorType,
           class AApproximationOperatorTypeLeft  = typename OperatorType::AOperatorType,
           class BFBTTypeRight                   = BFBTOperator< typename OperatorType::AOperatorType,
                                               AApproximationOperatorTypeRight,
                                               AApproximationOperatorTypeLeft,
                                               typename OperatorType::BTOperatorType,
                                               typename OperatorType::BOperatorType >,
           class BFBTTypeLeft                    = BFBTOperator< typename OperatorType::AOperatorType,
                                              AApproximationOperatorTypeRight,
                                              AApproximationOperatorTypeLeft,
                                              typename OperatorType::BTOperatorType,
                                              typename OperatorType::BOperatorType >,
           class BFBTSubOperatorSolverTypeRight  = hyteg::CGSolver< BFBTTypeRight >,
           class BFBTSubOperatorSolverTypeLeft   = hyteg::CGSolver< BFBTTypeLeft > >
class SchurBFBTSolver : public SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >
{
 public:
   using SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >::prefix_;
   typedef typename OperatorType::AOperatorType               AOperatorType;
   typedef typename OperatorType::BTOperatorType              BTOperatorType;
   typedef typename OperatorType::BOperatorType               BOperatorType;
   typedef typename OperatorType::srcType::PressureFunction_T PressureFunctionType;
   typedef BFBTOperator< typename OperatorType::AOperatorType,
                         AApproximationOperatorTypeRight,
                         AApproximationOperatorTypeLeft,
                         BTOperatorType,
                         BOperatorType >
       BFBTTypeInner;

   SchurBFBTSolver( walberla::Config::BlockHandle&                    parameters,
                    const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                    uint_t                                            minLevel,
                    uint_t                                            maxLevel,

                    const std::shared_ptr< OperatorType >&                                    saddleOp,
                    const std::shared_ptr< AApproximationOperatorTypeRight >&                 ApproximationOperatorRight,
                    const std::shared_ptr< AApproximationOperatorTypeLeft >&                  ApproximationOperatorLeft,
                    const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >& AApproximationSolverRight,
                    const std::shared_ptr< MantleConvection::ABlockSolver< AOperatorType > >& AApproximationSolverLeft,

                    hyteg::BoundaryCondition velocityBCx,
                    hyteg::BoundaryCondition velocityBCy,
                    hyteg::BoundaryCondition velocityBCz,

                    bool lowMemoryMode = false,

                    const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& BFBTPreconditionerRight =
                        std::make_shared< SchurIdentityPreconditioner< PressureFunctionType > >(),
                    const std::shared_ptr< BFBTTypeRight >&                                        BFBTOperatorRight = nullptr,
                    const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >& BFBTPreconditionerLeft =
                        std::make_shared< SchurIdentityPreconditioner< PressureFunctionType > >(),
                    const std::shared_ptr< BFBTTypeLeft >& BFBTOperatorLeft = nullptr,

                    std::string prefix = "" )
   : SchurSolver< SchurOperator< typename OperatorType::srcType::PressureFunction_T > >( prefix )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , lowMemoryMode_( lowMemoryMode )
   {
      auto substSolverARight_ =
          std::make_shared< hyteg::SubstitutePreconditioner< AApproximationOperatorTypeRight, AOperatorType > >(
              AApproximationSolverRight, saddleOp->getAPtr() );
      auto substSolverALeft_ =
          std::make_shared< hyteg::SubstitutePreconditioner< AApproximationOperatorTypeLeft, AOperatorType > >(
              AApproximationSolverLeft, saddleOp->getAPtr() );

      BFBTOperatorInner_ = std::make_shared< BFBTTypeInner >( storage,
                                                              minLevel,
                                                              maxLevel,
                                                              saddleOp->getAPtr(),
                                                              ApproximationOperatorRight,
                                                              ApproximationOperatorLeft,
                                                              saddleOp->getBTPtr(),
                                                              saddleOp->getBPtr(),
                                                              substSolverARight_,
                                                              substSolverALeft_,
                                                              velocityBCx,
                                                              velocityBCy,
                                                              velocityBCz,
                                                              lowMemoryMode_ );

      BFBTOperatorRight_ = BFBTOperatorRight;
      if constexpr ( std::is_same< BFBTTypeRight, BFBTTypeInner >::value )
      {
         if ( BFBTOperatorRight == nullptr )
         {
            BFBTOperatorRight_ = BFBTOperatorInner_;
         }
      }

      BFBTOperatorLeft_ = BFBTOperatorLeft;
      if constexpr ( std::is_same< BFBTTypeRight, BFBTTypeInner >::value )
      {
         if ( BFBTOperatorLeft == nullptr )
         {
            BFBTOperatorLeft_ = BFBTOperatorInner_;
         }
      }

      if ( BFBTOperatorRight_ == nullptr || BFBTOperatorLeft_ == nullptr )
      {
         WALBERLA_ABORT( "Right or left BFBT outer operator undefined." );
      }

      maxIterations_     = parameters.getParameter< uint_t >( prefix + std::string( "SchurBFBTMaxIterations" ) );
      relativeTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SchurBFBTRelativeTolerance" ) );
      absoluteTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "SchurBFBTAbsoluteTolerance" ) );
      printInfo_         = parameters.getParameter< bool >( prefix + std::string( "SchurBFBTPrintInfo" ) );

      // solvers
      auto substOperator_ = std::make_shared< SchurOperator< PressureFunctionType > >( storage, minLevel, maxLevel );

      // right solver
      auto substSolverRight_ =
          std::make_shared< hyteg::SubstitutePreconditioner< BFBTTypeRight, SchurOperator< PressureFunctionType > > >(
              BFBTPreconditionerRight, substOperator_ );

      if constexpr ( std::is_same< BFBTSubOperatorSolverTypeRight, hyteg::CGSolver< BFBTTypeRight > >::value )
      {
         OuterBFBTSolverRight_ = std::make_shared< BFBTSubOperatorSolverTypeRight >( storage,
                                                                                     minLevel,
                                                                                     maxLevel,
                                                                                     maxIterations_,
                                                                                     relativeTolerance_,
                                                                                     absoluteTolerance_,
                                                                                     substSolverRight_,
                                                                                     lowMemoryMode_ );
         OuterBFBTSolverRight_->setPrintInfo( printInfo_ );

         if ( prefix_ != "" )
         {
            OuterBFBTSolverRight_->setName( std::string( "CG Schur BFBT Right prefix " ) + prefix_ );
         }
         else
         {
            OuterBFBTSolverRight_->setName( "CG Schur BFBT Right" );
         }
      }
      else if constexpr ( std::is_same< BFBTSubOperatorSolverTypeRight, hyteg::FGMRESSolver< BFBTTypeRight > >::value )
      {
         OuterBFBTSolverRight_ = std::make_shared< BFBTSubOperatorSolverTypeRight >( storage,
                                                                                     minLevel,
                                                                                     maxLevel,
                                                                                     maxIterations_,
                                                                                     relativeTolerance_,
                                                                                     absoluteTolerance_,
                                                                                     substSolverRight_,
                                                                                     maxIterations_,
                                                                                     real_c( 0 ),
                                                                                     real_c( 0 ),
                                                                                     lowMemoryMode_ );
         OuterBFBTSolverRight_->setPrintInfo( printInfo_ );

         if ( prefix_ != "" )
         {
            OuterBFBTSolverRight_->setName( std::string( "FGMRES Schur BFBT Right prefix " ) + prefix_ );
         }
         else
         {
            OuterBFBTSolverRight_->setName( "FGMRES Schur BFBT Right" );
         }
      }
      else
      {
         WALBERLA_ABORT( "Unsupported BFBT Suboperator solver type." );
      }

      // left solver
      auto substSolverLeft_ =
          std::make_shared< hyteg::SubstitutePreconditioner< BFBTTypeLeft, SchurOperator< PressureFunctionType > > >(
              BFBTPreconditionerLeft, substOperator_ );

      if constexpr ( std::is_same< BFBTSubOperatorSolverTypeLeft, hyteg::CGSolver< BFBTTypeLeft > >::value )
      {
         OuterBFBTSolverLeft_ = std::make_shared< BFBTSubOperatorSolverTypeLeft >( storage,
                                                                                   minLevel,
                                                                                   maxLevel,
                                                                                   maxIterations_,
                                                                                   relativeTolerance_,
                                                                                   absoluteTolerance_,
                                                                                   substSolverLeft_,
                                                                                   lowMemoryMode_ );
         OuterBFBTSolverLeft_->setPrintInfo( printInfo_ );

         if ( prefix_ != "" )
         {
            OuterBFBTSolverLeft_->setName( std::string( "CG Schur BFBT Left prefix " ) + prefix_ );
         }
         else
         {
            OuterBFBTSolverLeft_->setName( "CG Schur BFBT Left" );
         }
      }
      else if constexpr ( std::is_same< BFBTSubOperatorSolverTypeLeft, hyteg::FGMRESSolver< BFBTTypeLeft > >::value )
      {
         OuterBFBTSolverLeft_ = std::make_shared< BFBTSubOperatorSolverTypeLeft >( storage,
                                                                                   minLevel,
                                                                                   maxLevel,
                                                                                   maxIterations_,
                                                                                   relativeTolerance_,
                                                                                   absoluteTolerance_,
                                                                                   substSolverRight_,
                                                                                   maxIterations_,
                                                                                   real_c( 0 ),
                                                                                   real_c( 0 ),
                                                                                   lowMemoryMode_ );
         OuterBFBTSolverLeft_->setPrintInfo( printInfo_ );

         if ( prefix_ != "" )
         {
            OuterBFBTSolverLeft_->setName( std::string( "FGMRES Schur BFBT Left prefix " ) + prefix_ );
         }
         else
         {
            OuterBFBTSolverLeft_->setName( "FGMRES Schur BFBT Left" );
         }
      }
      else
      {
         WALBERLA_ABORT( "Unsupported BFBT Suboperator solver type." );
      }

      // init functions
      if ( !lowMemoryMode_ )
      {
         temporary_         = std::make_shared< PressureFunctionType >( "bfbt temporary", storage, minLevel, maxLevel );
         temporarySolution_ = std::make_shared< PressureFunctionType >( "bfbt temporary solution", storage, minLevel, maxLevel );
      }
   }

   void solve( const SchurOperator< PressureFunctionType >& A,
               const PressureFunctionType&                  x,
               const PressureFunctionType&                  b,
               const walberla::uint_t                       level ) override
   {
      WALBERLA_UNUSED( A );

      std::shared_ptr< PressureFunctionType > temporaryApply;
      std::shared_ptr< PressureFunctionType > temporarySolutionApply;

      if ( !lowMemoryMode_ )
      {
         temporaryApply         = temporary_;
         temporarySolutionApply = temporarySolution_;
      }
      else
      {
         temporaryApply         = hyteg::getTemporaryFunction< PressureFunctionType >( storage_, minLevel_, maxLevel_ );
         temporarySolutionApply = hyteg::getTemporaryFunction< PressureFunctionType >( storage_, minLevel_, maxLevel_ );
      }

      temporaryApply->copyBoundaryConditionFromFunction( x );
      temporarySolutionApply->copyBoundaryConditionFromFunction( x );

      // Right solve
      if constexpr ( hyteg::SFINAE::has_setState< BFBTTypeRight >() )
      {
         BFBTOperatorRight_->setState( hyteg::BFBTState::RightBFBT );
      }
      temporaryApply->setToZero( level );
      OuterBFBTSolverRight_->solve( *BFBTOperatorRight_, *temporaryApply, b, level );

      // Inner solve
      if constexpr ( hyteg::SFINAE::has_setState< BFBTTypeInner >() )
      {
         BFBTOperatorInner_->setState( hyteg::BFBTState::InnerBFBT );
      }
      BFBTOperatorInner_->apply( *temporaryApply, *temporarySolutionApply, level, hyteg::All, hyteg::Replace );

      // Left solve
      if constexpr ( hyteg::SFINAE::has_setState< BFBTTypeLeft >() )
      {
         BFBTOperatorLeft_->setState( hyteg::BFBTState::LeftBFBT );
      }

      //x.setToZero( level );
      OuterBFBTSolverLeft_->solve( *BFBTOperatorLeft_, x, *temporarySolutionApply, level );
   };

   std::shared_ptr< BFBTSubOperatorSolverTypeRight > getSolverRight() { return OuterBFBTSolverRight_; }
   std::shared_ptr< BFBTSubOperatorSolverTypeLeft >  getSolverLeft() { return OuterBFBTSolverLeft_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#######################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "########## Schur BFBT Solver ##########"                                                                 << "\n";
      os << std::string( offset, ' ') << "#######################################"                                                                 << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "prefix_: "                     << prefix_                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "minLevel_: "                   << minLevel_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "maxLevel_: "                   << maxLevel_                  << "\n";            
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "maxIterations_: "              << maxIterations_             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "relativeTolerance_: "          << relativeTolerance_         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "absoluteTolerance_: "          << absoluteTolerance_         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "printInfo_: "                  << printInfo_                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 29 ) << std::left << "lowMemoryMode_: "              << lowMemoryMode_                    ;
      // clang-format on

      return os;
   }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< BFBTTypeInner > BFBTOperatorInner_;
   std::shared_ptr< BFBTTypeRight > BFBTOperatorRight_;
   std::shared_ptr< BFBTTypeLeft >  BFBTOperatorLeft_;

   std::shared_ptr< BFBTSubOperatorSolverTypeRight > OuterBFBTSolverRight_;
   std::shared_ptr< BFBTSubOperatorSolverTypeLeft >  OuterBFBTSolverLeft_;

   uint_t maxIterations_;
   real_t relativeTolerance_;
   real_t absoluteTolerance_;
   bool   printInfo_;
   bool   lowMemoryMode_;

   std::shared_ptr< PressureFunctionType > temporary_;
   std::shared_ptr< PressureFunctionType > temporarySolution_;
};

template < class MassOperatorType >
inline std::ostream& operator<<( std::ostream& os, const SchurBFBTSolver< MassOperatorType >& sdabs )
{
   return sdabs.print( os );
}

} // namespace MantleConvection