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
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

#include "ABlockIdentityPreconditioner.hpp"
#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class ABlockWeightedJacobiSolver : public ABlockSolver< OperatorType >
{
 public:
   using ABlockSolver< OperatorType >::prefix_;

   ABlockWeightedJacobiSolver( walberla::Config::BlockHandle&                    parameters,
                               const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
                               uint_t                                            minLevel,
                               uint_t                                            maxLevel,
                               std::string                                       prefix = "" )
   : ABlockSolver< OperatorType >( prefix )
   {
      iterations_ = parameters.getParameter< uint_t >( prefix + std::string( "ABlockWeightedJacobiIterations" ) );
      relaxParam_ = parameters.getParameter< real_t >( prefix + std::string( "ABlockWeightedJacobiRelaxParameter" ) );

      ABlockWeightedJacobiSolver_ =
          std::make_shared< hyteg::WeightedJacobiSmoother< OperatorType > >( storage, minLevel, maxLevel, relaxParam_ );
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      for ( uint_t i = 0; i < iterations_; i++ )
      {
         ABlockWeightedJacobiSolver_->solve( A, x, b, level );
      }
   };

   std::shared_ptr< hyteg::WeightedJacobiSmoother< OperatorType > > getSolver() { return ABlockWeightedJacobiSolver_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "######### ABlock Weighted Jacobi Solver #########"                                                       << "\n";
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "iterations_: "                << iterations_                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relaxParam_: "                << relaxParam_                        ;
      // clang-format on

      return os;
   }

 private:
   std::shared_ptr< hyteg::WeightedJacobiSmoother< OperatorType > > ABlockWeightedJacobiSolver_;
   uint_t                                                           iterations_;
   real_t                                                           relaxParam_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const ABlockWeightedJacobiSolver< OperatorType >& abwjs )
{
   return abwjs.print( os );
}

} // namespace MantleConvection