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
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

#include "SaddlePointSolver.hpp"

namespace MantleConvection {

template < class OperatorType >
class SaddlePointIdentityPreconditioner : public SaddlePointSolver< OperatorType >
{
 public:
   using SaddlePointSolver< OperatorType >::prefix_;

   SaddlePointIdentityPreconditioner()
   : SaddlePointSolver< OperatorType >()
   {
      SaddlePointIdentityPreconditioner_ = std::make_shared< hyteg::IdentityPreconditioner< OperatorType > >();
   }

   void solve( const OperatorType&                   A,
               const typename OperatorType::srcType& x,
               const typename OperatorType::dstType& b,
               const walberla::uint_t                level ) override
   {
      SaddlePointIdentityPreconditioner_->solve( A, x, b, level );
   };

   std::shared_ptr< hyteg::IdentityPreconditioner< OperatorType > > getSolver() { return SaddlePointIdentityPreconditioner_; }

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const override
   {
      // clang-format off
      os << std::string( offset, ' ') << "################################################"                    << "\n";
      os << std::string( offset, ' ') << "##### Saddle Point Identity Preconditioner #####"                    << "\n";
      os << std::string( offset, ' ') << "################################################"                    << "\n";
      os << std::string( offset, ' ') << "   " << "------No Parameters------"                                         ;
      // clang-format on

      return os;
   }

 private:
   std::shared_ptr< hyteg::IdentityPreconditioner< OperatorType > > SaddlePointIdentityPreconditioner_;
};

template < class OperatorType >
inline std::ostream& operator<<( std::ostream& os, const SaddlePointIdentityPreconditioner< OperatorType >& spip )
{
   return spip.print( os );
}

} // namespace MantleConvection