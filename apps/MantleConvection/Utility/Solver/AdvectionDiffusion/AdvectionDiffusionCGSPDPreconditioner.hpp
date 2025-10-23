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

#include "../../LHS/AdvectionDiffusionSPDSubOperator.hpp"
#include "AdvectionDiffusionIdentityPreconditioner.hpp"
#include "AdvectionDiffusionSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class AdvectionDiffusionCGSPDPreconditioner : public AdvectionDiffusionSolver< OperatorType >
{
 public:
   using AdvectionDiffusionSolver< OperatorType >::prefix_;

   AdvectionDiffusionCGSPDPreconditioner(
       walberla::Config::BlockHandle&                    parameters,
       const std::shared_ptr< hyteg::PrimitiveStorage >& storage,
       uint_t                                            minLevel,
       uint_t                                            maxLevel,
       const std::shared_ptr< OperatorType >&            AdvectionDiffusionOperator,
       const std::shared_ptr< AdvectionDiffusionSolver< AdvectionDiffusionSPDSubOperator< OperatorType > > >&
           AdvectionDiffusionPreconditioner =
               std::make_shared< AdvectionDiffusionIdentityPreconditioner< AdvectionDiffusionSPDSubOperator< OperatorType > > >(),
       std::string prefix = "" )
   : AdvectionDiffusionSolver< OperatorType >( prefix )
   , AdvectionDiffusionPreconditioner_( AdvectionDiffusionPreconditioner )
   , AdvectionDiffusionOperator_( AdvectionDiffusionOperator )
   {
      maxIterations_ =
          parameters.getParameter< uint_t >( prefix + std::string( "AdvectionDiffusionCGSPDPreconditionerMaxIterations" ) );
      relativeTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "AdvectionDiffusionCGSPDPreconditionerRelativeTolerance" ) );
      absoluteTolerance_ =
          parameters.getParameter< real_t >( prefix + std::string( "AdvectionDiffusionCGSPDPreconditionerAbsoluteTolerance" ) );
      printInfo_ = parameters.getParameter< bool >( prefix + std::string( "AdvectionDiffusionCGSPDPreconditionerPrintInfo" ) );

      AdvectionDiffusionSPDSubOperator_ = std::make_shared< AdvectionDiffusionSPDSubOperator< OperatorType > >(
          storage, minLevel, maxLevel, AdvectionDiffusionOperator_ );

      AdvectionDiffusionCGSPDPreconditioner_ =
          std::make_shared< hyteg::CGSolver< AdvectionDiffusionSPDSubOperator< OperatorType > > >(
              storage,
              minLevel,
              maxLevel,
              maxIterations_,
              relativeTolerance_,
              absoluteTolerance_,
              AdvectionDiffusionPreconditioner_ );
      AdvectionDiffusionCGSPDPreconditioner_->setPrintInfo( printInfo_ );

      if ( prefix_ != "" )
      {
         AdvectionDiffusionCGSPDPreconditioner_->setName( std::string( "CG Advection Diffusion CG SPD Preconditioner prefix " ) +
                                                          prefix_ );
      }
      else
      {
         AdvectionDiffusionCGSPDPreconditioner_->setName( "CG Advection Diffusion CG SPD Preconditioner" );
      }
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      WALBERLA_UNUSED( A );
      AdvectionDiffusionCGSPDPreconditioner_->solve( *AdvectionDiffusionSPDSubOperator_, x, b, level );
   };

   std::shared_ptr< hyteg::CGSolver< AdvectionDiffusionSPDSubOperator< OperatorType > > > getSolver()
   {
      return AdvectionDiffusionCGSPDPreconditioner_;
   }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "###################################################"                                                     << "\n";
      os << std::string( offset, ' ') << "#### Advection Diffusion CG SPD Preconditioner ####"                                                     << "\n";
      os << std::string( offset, ' ') << "###################################################"                                                     << "\n";
      os << std::string( offset, ' ') << "   " << "------Parameters------"                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "prefix_: "                    << prefix_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "maxIterations_: "             << maxIterations_              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "relativeTolerance_: "         << relativeTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "absoluteTolerance_: "         << absoluteTolerance_          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 21 ) << std::left << "printInfo_: "                 << printInfo_                  << "\n";
      // clang-format on

      AdvectionDiffusionPreconditioner_->print( os, offset + 3 );

      return os;
   }

 private:
   std::shared_ptr< hyteg::CGSolver< AdvectionDiffusionSPDSubOperator< OperatorType > > > AdvectionDiffusionCGSPDPreconditioner_;
   std::shared_ptr< AdvectionDiffusionSolver< AdvectionDiffusionSPDSubOperator< OperatorType > > >
                                                                       AdvectionDiffusionPreconditioner_;
   std::shared_ptr< OperatorType >                                     AdvectionDiffusionOperator_;
   std::shared_ptr< AdvectionDiffusionSPDSubOperator< OperatorType > > AdvectionDiffusionSPDSubOperator_;

   uint_t maxIterations_;
   real_t relativeTolerance_;
   real_t absoluteTolerance_;
   bool   printInfo_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const AdvectionDiffusionCGSPDPreconditioner< OperatorType >& abcgols )
{
   return abcgols.print( os );
}

} // namespace MantleConvection