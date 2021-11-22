/*
 * Copyright (c) 2020 Marcus Mohr.
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

#include <typeinfo>

#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/indexing/CouplingCount.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {
namespace indexing {

/// Return the number of couplings between all DoFs for the given operator
///
/// The function determines the couplings between all DoFs (including those that are fixed by Dirichlet
/// boundary conditions) for the given operator. The number equals the number of non-zero entries, if
/// the matrix corresponding to the operator is assembled.
///
/// \note: The determination of the couplings is based on the structure of the operator and the
///        supports of the ansatz function in its associated spaced. This is to say the function
///        does not check if coupling values equal zero numerically.
///
/// \note: The function can directly handle all scalar P1 and P2 operators, as here the couplings
///        are determined solely from the supports. In the case of composite (block) operators
///        this must currently be done manually and _all_ operators must be added to the corresponding
///        if-branch.
template < typename opType >
uint_t getNumberOfGlobalDoFCouplings( const opType& oper, uint_t level )
{
   uint_t nCouplings = 0;

   typedef Operator< P1Function< real_t >, P1Function< real_t > >           P1ScalarOp;
   typedef Operator< P2Function< real_t >, P2Function< real_t > >           P2ScalarOp;
   typedef Operator< EdgeDoFFunction< real_t >, EdgeDoFFunction< real_t > > EdgeDoFScalarOp;

   // Scalar P1 operators
   if ( dynamic_cast< const P1ScalarOp* >( &oper ) )
   {
      nCouplings = countLocalDoFCouplings< P1FunctionTag, P1FunctionTag >( oper.getStorage(), level );
   }

   // Scalar P2 operators
   else if ( dynamic_cast< const P2ScalarOp* >( &oper ) )
   {
      nCouplings = countLocalDoFCouplings< P2FunctionTag, P2FunctionTag >( oper.getStorage(), level );
   }

   // Scalar EdgeDoF operators
   else if ( dynamic_cast< const EdgeDoFScalarOp* >( &oper ) )
   {
      nCouplings = countLocalDoFCouplings< EdgeDoFFunctionTag, EdgeDoFFunctionTag >( oper.getStorage(), level );
   }

   // P2-P1 Taylor-Hood pure Stokes variants
   else if ( dynamic_cast< const P2P1TaylorHoodStokesOperator* >( &oper ) ||
             dynamic_cast< const P2P1ElementwiseBlendingStokesOperator* >( &oper ) )
   {
      uint_t dim = oper.getStorage()->hasGlobalCells() ? 3 : 2;
      nCouplings = dim * countLocalDoFCouplings< P2FunctionTag, P2FunctionTag >( oper.getStorage(), level ) +
                   2 * countLocalDoFCouplings< P2VectorFunctionTag, P1FunctionTag >( oper.getStorage(), level );
   }
   else
   {
      WALBERLA_ABORT( "Operator resolution failure in indexing::getNumberOfGlobalDoFCouplings()" );
   }

   // MPI communication
   walberla::mpi::allReduceInplace( nCouplings, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );

   return nCouplings;
}

} // namespace indexing
} // namespace hyteg
