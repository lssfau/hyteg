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
#include "hyteg/solvers/Solver.hpp"

#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
class ABlockFullMultigridSolver : public ABlockSolver< OperatorType >
{
 public:
   using ABlockSolver< OperatorType >::prefix_;

   ABlockFullMultigridSolver(
       walberla::Config::BlockHandle&                                           parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                        storage,
       uint_t                                                                   minLevel,
       uint_t                                                                   maxLevel,
       const std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > >& ABlockMGSolver,
       const std::shared_ptr< RestrictionOperatorType >&                        ABlockRestriction,
       const std::shared_ptr< ProlongationOperatorType >&                       ABlockProlongation,
       std::function< void( const typename OperatorType::srcType&, const typename OperatorType::dstType&, uint_t ) >
           preSolveCallback =
               []( const typename OperatorType::srcType& x, const typename OperatorType::dstType& b, uint_t level ) {
                  WALBERLA_UNUSED( x );
                  WALBERLA_UNUSED( level );
               },
       hyteg::DoFType flag   = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
       std::string    prefix = "" )
   : ABlockSolver< OperatorType >( prefix )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , ABlockMGSolver_( ABlockMGSolver )
   , ABlockRestriction_( ABlockRestriction )
   , ABlockProlongation_( ABlockProlongation )
   , preSolveCallback_( preSolveCallback )
   , flag_( flag )
   {
      cyclesPerLevel_ = parameters.getParameter< uint_t >( prefix + std::string( "ABlockFullMultigridCyclesPerLevel" ) );
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      for ( uint_t currentLevel = level; currentLevel > minLevel_; currentLevel-- )
      {
         ABlockRestriction_->restrict( b, level, hyteg::All );
      }

      for ( uint_t currentLevel = minLevel_; currentLevel <= level; currentLevel++ )
      {
         preSolveCallback_( x, b, currentLevel );

         for ( uint_t k = 0; k < cyclesPerLevel_; k++ )
         {
            ABlockMGSolver_->solve( A, x, b, currentLevel );
         }

         if ( currentLevel < level )
         {
            ABlockProlongation_->prolongate( x, currentLevel, flag_ );
         }
      }
   };

   std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > > getSolver()
   {
      return ABlockMGSolver_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "######################################################"                                         << "\n";
      os << std::string( offset, ' ') << "############ ABlock Full Multigrid Solver ############"                                         << "\n";
      os << std::string( offset, ' ') << "######################################################"                                         << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "prefix_: "                    << prefix_            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "cyclesPerLevel_: "            << cyclesPerLevel_    << "\n";
      // clang-format on

      ABlockMGSolver_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > > ABlockMGSolver_;
   std::shared_ptr< RestrictionOperatorType >                        ABlockRestriction_;
   std::shared_ptr< ProlongationOperatorType >                       ABlockProlongation_;

   std::function< void( const typename OperatorType::srcType&, const typename OperatorType::dstType&, uint_t ) >
       preSolveCallback_;

   uint_t cyclesPerLevel_;

   hyteg::DoFType flag_;
};

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
inline std::ostream&
    operator<<( std::ostream&                                                                                       os,
                const ABlockFullMultigridSolver< OperatorType, RestrictionOperatorType, ProlongationOperatorType >& abfmgs )
{
   return abfmgs.print( os );
}

} // namespace MantleConvection