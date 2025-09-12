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

#include "ABlockIdentityPreconditioner.hpp"
#include "ABlockSolver.hpp"

namespace MantleConvection {

template < class OperatorType, class MassOperatorType >
class ABlockCGVectorMassSolver : public ABlockSolver< OperatorType >
{
 public:
   using ABlockSolver< OperatorType >::prefix_;

   ABlockCGVectorMassSolver( walberla::Config::BlockHandle&                             parameters,
                             const std::shared_ptr< hyteg::PrimitiveStorage >&          storage,
                             uint_t                                                     minLevel,
                             uint_t                                                     maxLevel,
                             const std::shared_ptr< MassOperatorType >&                 MassOperator,
                             const std::shared_ptr< ABlockSolver< MassOperatorType > >& MassPreconditioner =
                                 std::make_shared< ABlockIdentityPreconditioner< MassOperatorType > >(),
                             bool        lowMemoryMode = false,
                             std::string prefix        = "" )
   : ABlockSolver< OperatorType >( prefix )
   , MassOperator_( MassOperator )
   , MassPreconditioner_( MassPreconditioner )
   , lowMemoryMode_( lowMemoryMode )
   {
      maxIterations_     = parameters.getParameter< uint_t >( prefix + std::string( "ABlockCGVectorMassMaxIterations" ) );
      relativeTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "ABlockCGVectorMassRelativeTolerance" ) );
      absoluteTolerance_ = parameters.getParameter< real_t >( prefix + std::string( "ABlockCGVectorMassAbsoluteTolerance" ) );
      printInfo_         = parameters.getParameter< bool >( prefix + std::string( "ABlockCGVectorMassPrintInfo" ) );

      MassSolver_ = std::make_shared< hyteg::CGSolver< MassOperatorType > >( storage,
                                                                             minLevel,
                                                                             maxLevel,
                                                                             maxIterations_,
                                                                             relativeTolerance_,
                                                                             absoluteTolerance_,
                                                                             MassPreconditioner_,                                                                             
                                                                             lowMemoryMode_ );
      MassSolver_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         MassSolver_->setName( std::string( "CG ABlock Vector Mass prefix " ) + prefix_ );
      }
      else
      {
         MassSolver_->setName( "CG ABlock Vector Mass" );
      }
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      WALBERLA_UNUSED( A );
      // WALBERLA_LOG_INFO_ON_ROOT( "ABlock CG Vector MassSolver solving..." );
      MassSolver_->solve( *MassOperator_, x, b, level );
   };

   std::shared_ptr< hyteg::CGSolver< MassOperatorType > > getSolver()
   {
      return MassSolver_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "########## ABlock CG Vector MassSolver ##########"                                                       << "\n";
      os << std::string( offset, ' ') << "#################################################"                                                       << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxIterations_: "             << maxIterations_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relativeTolerance_: "         << relativeTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "printInfo_: "                 << printInfo_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "lowMemoryMode_: "             << lowMemoryMode_              << "\n";
      // clang-format on

      return os;
   }

 private:
   std::shared_ptr< hyteg::CGSolver< MassOperatorType > > MassSolver_;
   std::shared_ptr< MassOperatorType >                    MassOperator_;
   std::shared_ptr< ABlockSolver< MassOperatorType > >    MassPreconditioner_;
   uint_t                                                 maxIterations_;
   real_t                                                 relativeTolerance_;
   real_t                                                 absoluteTolerance_;
   bool                                                   printInfo_;

   bool lowMemoryMode_;
};

template < class OperatorType, class MassOperatorType >
inline std::ostream& operator<<( std::ostream& os, const ABlockCGVectorMassSolver< OperatorType, MassOperatorType >& abcgvms )
{
   return abcgvms.print( os );
}

} // namespace MantleConvection