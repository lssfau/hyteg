/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#include "core/timing/Timer.h"

#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/solvers/Solver.hpp"

namespace hyteg {

/// \brief Solver wrapper that implements agglomeration.
///
/// Agglomeration can be effective in massively parallel settings if
/// a solver does not scale well and should be called by only a subset
/// of processes.
///
/// For example this can make sense in a multigrid iteration as the
/// coarse grid solver usually solves a much smaller problem.
///
/// Therefore, the problem can be re-distributed to a subset of processes
/// before the coarse grid solver is called and the solution is then
/// distributed back to the original processes.
///
/// This wrapper can be used to wrap the coarse grid solver so that
/// the agglomeration process can be implemented conveniently.
///
/// The following steps must be performed to perform agglomeration (exact code might change):
///
///     1. a copy of the PrimitiveStorage must be created
///          auto agglomerationStorage = storage->createCopy();
///
///     2. the coarse grid solver must be allocated using the copied storage
///          auto coarseGridSolver = std::make_shared< CGSolver< ... > >( agglomerationStorage, ... );
///
///     3. the coarse grid solver must be wrapped with this wrapper
///          auto wrappedCoarseGridSolver = std::make_shared< AgglomerationWrapper< ... > >( coarseGridSolver, agglomerationStorage ... );
///
///     4. the wrapper is passed as the coarse grid solver to the multigrid solver
///          auto mgSolver = std::make_shared< MGSolver< ... > >( storage, wrappedCoarseGridSolver, ... );
///
/// Before every coarse grid solve, the wrapper then distributes the corresponding functions to
/// a subset of processes, calls the solver, and then distributes the solution back to the original
/// storage.
///
template < typename OperatorType >
class AgglomerationWrapper : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   AgglomerationWrapper( const std::shared_ptr< Solver< OperatorType > >& solver,
                         const std::shared_ptr< PrimitiveStorage >&       agglomerationStorage,
                         const uint_t&                                    level,
                         const uint_t&                                    numberOfAgglomerationProcesses )
   : solver_( solver )
   , agglomerationStorage_( agglomerationStorage )
   , level_( level )
   , numberOfAgglomerationProcesses_( numberOfAgglomerationProcesses )
   , A_agglomeration( agglomerationStorage, level, level )
   , b_agglomeration( "b_agglomeration", agglomerationStorage, level, level )
   , x_agglomeration( "x_agglomeration", agglomerationStorage, level, level )
   {}

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_CHECK_EQUAL( level, level_, "Agglomeration was prepared for level " << level_ );

      // asserting that the Primitive distribution is equal here
      // copying the distribution would be more costly but safer
      WALBERLA_CHECK( x.getStorage()->getPrimitiveIDs() == agglomerationStorage_->getPrimitiveIDs() );

      b_agglomeration.copyFrom( b, level );
      x_agglomeration.copyFrom( x, level );

      loadbalancing::distributed::roundRobin( *agglomerationStorage_, numberOfAgglomerationProcesses_ );

      solver_->solve( A_agglomeration, x_agglomeration, b_agglomeration, level );

      loadbalancing::distributed::copyDistribution( *x.getStorage(), *agglomerationStorage_ );

      x.copyFrom( x_agglomeration, level );
   }

 private:
   std::shared_ptr< Solver< OperatorType > > solver_;
   std::shared_ptr< PrimitiveStorage >       agglomerationStorage_;
   uint_t                                    level_;
   uint_t                                    numberOfAgglomerationProcesses_;

   OperatorType A_agglomeration;
   FunctionType b_agglomeration;
   FunctionType x_agglomeration;
};

} // namespace hyteg