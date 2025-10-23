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

#include "SchurSolver.hpp"

namespace MantleConvection {

template < class PressureFunctionType, class RestrictionOperatorType, class ProlongationOperatorType >
class SchurFullMultigridSolver : public SchurSolver< SchurOperator< PressureFunctionType > >
{
 public:
   using SchurSolver< SchurOperator< PressureFunctionType > >::prefix_;

   SchurFullMultigridSolver(
       walberla::Config::BlockHandle&                                                            parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >&                                         storage,
       uint_t                                                                                    minLevel,
       uint_t                                                                                    maxLevel,
       const std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > >&            SchurMGSolver,
       const std::shared_ptr< RestrictionOperatorType >&                                         SchurRestriction,
       const std::shared_ptr< ProlongationOperatorType >&                                        SchurProlongation,
       std::function< void( const PressureFunctionType&, const PressureFunctionType&, uint_t ) > preSolveCallback =
           []( const PressureFunctionType& x, const PressureFunctionType& b, uint_t level ) {
              WALBERLA_UNUSED( x );
              WALBERLA_UNUSED( level );
           },
       hyteg::DoFType flag   = hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary,
       std::string    prefix = "" )
   : SchurSolver< SchurOperator< PressureFunctionType > >( prefix )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , SchurMGSolver_( SchurMGSolver )
   , SchurRestriction_( SchurRestriction )
   , SchurProlongation_( SchurProlongation )
   , preSolveCallback_( preSolveCallback )
   , flag_( flag )
   {
      cyclesPerLevel_ = parameters.getParameter< uint_t >( prefix + std::string( "SchurFullMultigridCyclesPerLevel" ) );
   }

   void solve( const SchurOperator< PressureFunctionType >& A,
               const PressureFunctionType&                  x,
               const PressureFunctionType&                  b,
               const walberla::uint_t                       level ) override
   {
      for ( uint_t currentLevel = level; currentLevel > minLevel_; currentLevel-- )
      {
         SchurRestriction_->restrict( b, level, hyteg::All );
      }

      for ( uint_t currentLevel = minLevel_; currentLevel <= level; currentLevel++ )
      {
         preSolveCallback_( x, b, currentLevel );

         for ( uint_t k = 0; k < cyclesPerLevel_; k++ )
         {
            SchurMGSolver_->solve( A, x, b, currentLevel );
         }

         if ( currentLevel < level )
         {
            SchurProlongation_->prolongate( x, currentLevel, flag_ );
         }
      }
   };

   std::shared_ptr< MantleConvection::ABlockSolver< OperatorType > > getSolver() { return SchurMGSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#####################################################"                                          << "\n";
      os << std::string( offset, ' ') << "############ Schur Full Multigrid Solver ############"                                          << "\n";
      os << std::string( offset, ' ') << "#####################################################"                                          << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "prefix_: "                    << prefix_            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 27 ) << std::left << "cyclesPerLevel_: "            << cyclesPerLevel_    << "\n";
      // clang-format on

      SchurMGSolver_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< SchurSolver< SchurOperator< PressureFunctionType > > > SchurMGSolver_;
   std::shared_ptr< RestrictionOperatorType >                              SchurRestriction_;
   std::shared_ptr< ProlongationOperatorType >                             SchurProlongation_;

   std::function< void( const PressureFunctionType&, const PressureFunctionType&, uint_t ) > preSolveCallback_;

   uint_t cyclesPerLevel_;

   hyteg::DoFType flag_;
};

template < class OperatorType, class RestrictionOperatorType, class ProlongationOperatorType >
inline std::ostream&
    operator<<( std::ostream&                                                                                      os,
                const SchurFullMultigridSolver< OperatorType, RestrictionOperatorType, ProlongationOperatorType >& sfmgs )
{
   return sfmgs.print( os );
}

} // namespace MantleConvection