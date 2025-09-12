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

#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/controlflow/AgglomerationWrapper.hpp"

#include "ABlockIdentityPreconditioner.hpp"
#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class ABlockAgglomeratedCGCoarseGridSolver : public ABlockSolver< OperatorType >
{
 public:
   using ABlockSolver< OperatorType >::prefix_;

   ABlockAgglomeratedCGCoarseGridSolver( walberla::Config::BlockHandle&                         parameters,
                                         const std::shared_ptr< hyteg::PrimitiveStorage >&      storage,
                                         uint_t                                                 minLevel,
                                         uint_t                                                 maxLevel,
                                         bool                                                   lowMemoryMode = false,
                                         const std::shared_ptr< ABlockSolver< OperatorType > >& ABlockPreconditioner =
                                             std::make_shared< ABlockIdentityPreconditioner< OperatorType > >(),
                                         std::string prefix = "" )
   : ABlockSolver< OperatorType >( prefix )
   , ABlockPreconditioner_( ABlockPreconditioner )
   , auxiliaryUpdateFunction_( []() {} )
   , lowMemoryMode_( lowMemoryMode )
   {
      maxIterations_ = parameters.getParameter< uint_t >( prefix + std::string( "ABlockAgglomeratedCGCoarseGridMaxIterations" ) );
      relativeTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "ABlockAgglomeratedCGCoarseGridRelativeTolerance" ) );
      absoluteTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "ABlockAgglomeratedCGCoarseGridAbsoluteTolerance" ) );
      printInfo_ = parameters.getParameter< bool >( prefix + std::string( "ABlockAgglomeratedCGCoarseGridPrintInfo" ) );

      if ( storage->hasGlobalCells() )
      {
         nMacroElementsGlobal_ = storage->getNumberOfGlobalCells();
      }
      else
      {
         nMacroElementsGlobal_ = storage->getNumberOfGlobalFaces();
      }

      uint_t paramRank = parameters.getParameter< uint_t >( prefix + std::string( "ABlockAgglomeratedCGCoarseGridRanks" ) );

      maxRank_ = std::min( nMacroElementsGlobal_, std::max( walberla::numeric_cast< uint_t >( 1 ), paramRank ) );

      agglomerationStorage_ = storage->createCopy();

      migrationInfoToAgglomerationStorage_ = loadbalancing::distributed::roundRobin( *agglomerationStorage_, 0, maxRank_ - 1 );
      WALBERLA_MPI_BARRIER();
      migrationInfoToOriginalStorage_ =
          loadbalancing::distributed::reverseDistributionDry( migrationInfoToAgglomerationStorage_ );
      WALBERLA_MPI_BARRIER();

      x_agglomeration_ =
          std::make_shared< typename OperatorType::srcType >( "xAgglomeration", agglomerationStorage_, minLevel, minLevel );
      b_agglomeration_ =
          std::make_shared< typename OperatorType::dstType >( "bAgglomeration", agglomerationStorage_, minLevel, minLevel );

      ABlockAgglomeratedCGCoarseGridSolver_ = std::make_shared< hyteg::CGSolver< OperatorType > >( agglomerationStorage_,
                                                                                                   minLevel,
                                                                                                   maxLevel,
                                                                                                   maxIterations_,
                                                                                                   relativeTolerance_,
                                                                                                   absoluteTolerance_,
                                                                                                   ABlockPreconditioner_,
                                                                                                   lowMemoryMode_ );
      ABlockAgglomeratedCGCoarseGridSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         ABlockAgglomeratedCGCoarseGridSolver_->setName( std::string( "CG ABlock Coarse Grid prefix " ) + prefix_ );
      }
      else
      {
         ABlockAgglomeratedCGCoarseGridSolver_->setName( "CG ABlock Coarse Grid" );
      }
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      WALBERLA_UNUSED( A );

      copyFunctionToAgglomerationStorage( b, *b_agglomeration_, level );
      copyFunctionToAgglomerationStorage( x, *x_agglomeration_, level );

      x_agglomeration_->copyBoundaryConditionFromFunction( x );

      auxiliaryUpdateFunction_();

      ABlockAgglomeratedCGCoarseGridSolver_->solve( *allgomeratedOperator_, *x_agglomeration_, *b_agglomeration_, level );

      copyFunctionToOriginalStorage( x, *x_agglomeration_, level );
   };

   std::shared_ptr< hyteg::CGSolver< OperatorType > > getSolver() { return ABlockAgglomeratedCGCoarseGridSolver_; }
   std::shared_ptr< hyteg::PrimitiveStorage >         getAgglomerationStorage() { return agglomerationStorage_; }

   template < class FunctionType >
   void copyFunctionToAgglomerationStorage( const FunctionType& f, const FunctionType& fAgglomerated, uint_t level )
   {
      fAgglomerated.copyFrom( f, level, migrationInfoToOriginalStorage_.getMap(), migrationInfoToAgglomerationStorage_.getMap() );
   }

   template < class FunctionType >
   void copyFunctionToOriginalStorage( const FunctionType& f, const FunctionType& fAgglomerated, uint_t level )
   {
      f.copyFrom( fAgglomerated, level, migrationInfoToAgglomerationStorage_.getMap(), migrationInfoToOriginalStorage_.getMap() );
   }

   void setAgglomerationOperator( const std::shared_ptr< OperatorType >& allgomeratedOperator )
   {
      allgomeratedOperator_ = allgomeratedOperator;
   }

   void setAuxiliaryUpdateFunction( std::function< void() > auxiliaryUpdateFunction )
   {
      auxiliaryUpdateFunction_ = auxiliaryUpdateFunction;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#######################################################"                                                 << "\n";
      os << std::string( offset, ' ') << "###### ABlock Agglomerated CG Coarse Grid Solver ######"                                                 << "\n";
      os << std::string( offset, ' ') << "#######################################################"                                                 << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxIterations_: "             << maxIterations_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relativeTolerance_: "         << relativeTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "printInfo_: "                 << printInfo_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "minRank_: "                   << minRank_                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxRank_: "                   << maxRank_                    << "\n";

      // clang-format on

      ABlockPreconditioner_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::CGSolver< OperatorType > > ABlockAgglomeratedCGCoarseGridSolver_;
   std::shared_ptr< ABlockSolver< OperatorType > >    ABlockPreconditioner_;
   uint_t                                             maxIterations_;
   real_t                                             relativeTolerance_;
   real_t                                             absoluteTolerance_;
   bool                                               printInfo_;

   std::shared_ptr< hyteg::PrimitiveStorage > agglomerationStorage_;
   hyteg::MigrationInfo                       migrationInfoToAgglomerationStorage_;
   hyteg::MigrationInfo                       migrationInfoToOriginalStorage_;

   std::shared_ptr< OperatorType >                   allgomeratedOperator_;
   std::shared_ptr< typename OperatorType::dstType > b_agglomeration_;
   std::shared_ptr< typename OperatorType::srcType > x_agglomeration_;

   std::function< void() > auxiliaryUpdateFunction_;

   uint_t nMacroElementsGlobal_;
   uint_t minRank_;
   uint_t maxRank_;

   bool lowMemoryMode_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const ABlockAgglomeratedCGCoarseGridSolver< OperatorType >& abcgs )
{
   return abcgs.print( os );
}

} // namespace MantleConvection