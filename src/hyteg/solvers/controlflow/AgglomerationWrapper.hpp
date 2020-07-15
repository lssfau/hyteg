/*
 * Copyright (c) 2017-2020 Nils Kohl.
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
/// Therefore, the problem is re-distributed to a subset of processes
/// before the coarse grid solver is called and the solution is then
/// distributed back to the original processes.
///
/// This wrapper can be used to wrap the coarse grid solver so that
/// the agglomeration process can be implemented conveniently.
///
/// The following steps must be performed to perform agglomeration (exact code might change):
///
/// 1. Create an instance of the agglomeration wrapper using the initial, globally distributed storage:
///
///        auto agglomerationWrapper = std::make_shared< AgglomerationWrapper >( storage, level, ... );
///
/// 2. Choose an agglomeration strategy;
///
///        agglomerationWrapper->setStrategyContinuousProcesses( 0, 23 );
///
/// 3. Internally, a copy of the storage is created and distributed to a subset of processes. Use this storage to assemble
///    the solver, e.g.:
///
///        auto coarseGridSolver = std::make_shared< MinResSolver >( agglomerationWrapper->getAgglomerationStorage(), ... );
///        agglomerationWrapper->setSolver( coarseGridSolver );
///
///    Note: If you require the operator of the coarse grid solver (e.g. to perform the factorization in an extra step),
///          remember to use the operator of the agglomeration wrapper via the getter, e.g.:
///
///
///
/// 4. Now the AgglomerationWrapper instance can be used as a coarse grid solver in massively parallel settings.
///
/// Notes:
///     - do not re-partition the agglomeration storage manually, the partitioning is done by the wrapper
///     - do not re-partition the original storage after creating this wrapper
///     - currently the operator is constructed internally, a setter can be implemented if necessary
///
template < typename OperatorType >
class AgglomerationWrapper : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   /// \brief Constructs the agglomeration wrapper.
   ///
   /// \param originalStorage the original PrimitiveStorage instance
   /// \param level the refinement level that is subject to agglomeration (in MG settings usually the coarsest level)
   /// \param solveOnEmptyProcesses if false, the solve() call is only performed on processes that own primitives,
   ///                              this might speed up certain solvers and save resources on empty processes, however,
   ///                              some solvers might not work in this case and deadlock, use with care
   AgglomerationWrapper( const std::shared_ptr< PrimitiveStorage >& originalStorage,
                         const uint_t&                              level,
                         const bool&                                solveOnEmptyProcesses = true )
   : originalStorage_( originalStorage )
   , level_( level )
   , solveOnEmptyProcesses_( solveOnEmptyProcesses )
   , isStrategySet_( false )
   , isAgglomerationProcess_( false )
   {}

   void setStrategyContinuousProcesses( const uint_t& minRank, const uint_t& maxRank )
   {
      WALBERLA_CHECK( !isStrategySet_ );

      agglomerationStorage_                = originalStorage_->createCopy();
      migrationInfoToAgglomerationStorage_ = loadbalancing::distributed::roundRobin( *agglomerationStorage_, minRank, maxRank );
      isAgglomerationProcess_              = uint_c( walberla::mpi::MPIManager::instance()->rank() ) >= minRank &&
                                uint_c( walberla::mpi::MPIManager::instance()->rank() ) <= maxRank;
      finalizeAgglomerationStrategy( migrationInfoToAgglomerationStorage_ );
   }

   void setStrategyEveryNthProcess( const uint_t& interval, const uint_t& numProcesses )
   {
      WALBERLA_CHECK( !isStrategySet_ );

      agglomerationStorage_ = originalStorage_->createCopy();
      migrationInfoToAgglomerationStorage_ =
          loadbalancing::distributed::roundRobinInterval( *agglomerationStorage_, interval, numProcesses );
      isAgglomerationProcess_ = uint_c( walberla::mpi::MPIManager::instance()->rank() ) % interval == 0 &&
                                uint_c( walberla::mpi::MPIManager::instance()->rank() ) / interval < numProcesses;
      finalizeAgglomerationStrategy( migrationInfoToAgglomerationStorage_ );
   }

   std::shared_ptr< PrimitiveStorage > getAgglomerationStorage() const { return agglomerationStorage_; }

   std::shared_ptr< OperatorType > getAgglomerationOperator() const { return A_agglomeration_; }

   bool isAgglomerationProcess() const
   {
      WALBERLA_CHECK( isStrategySet_ );
      return isAgglomerationProcess_;
   }

   void setSolver( const std::shared_ptr< Solver< OperatorType > >& solver ) { solver_ = solver; }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      WALBERLA_CHECK( isStrategySet_, "An agglomeration strategy must be set explicitly before solving." )
      WALBERLA_CHECK_EQUAL( level, level_, "Agglomeration was prepared for level " << level_ );
      WALBERLA_CHECK_NOT_NULLPTR( solver_.get(), "Solver for AgglomerationWrapper was not set." )

      const bool emptyProcess = agglomerationStorage_->getNumberOfLocalPrimitives() == 0;

      WALBERLA_MPI_BARRIER();

      b_agglomeration_->copyFrom(
          b, level, migrationInfoToOriginalStorage_.getMap(), migrationInfoToAgglomerationStorage_.getMap() );
      x_agglomeration_->copyFrom(
          x, level, migrationInfoToOriginalStorage_.getMap(), migrationInfoToAgglomerationStorage_.getMap() );

      if ( solveOnEmptyProcesses_ || !emptyProcess )
      {
         solver_->solve( *A_agglomeration_, *x_agglomeration_, *b_agglomeration_, level );
      }

      x.copyFrom(
          *x_agglomeration_, level, migrationInfoToAgglomerationStorage_.getMap(), migrationInfoToOriginalStorage_.getMap() );
   }

 private:
   void finalizeAgglomerationStrategy( const MigrationInfo& originalMigrationInfo )
   {
      migrationInfoToOriginalStorage_ = loadbalancing::distributed::reverseDistributionDry( originalMigrationInfo );

      A_agglomeration_ = std::make_shared< OperatorType >( agglomerationStorage_, level_, level_ );
      x_agglomeration_ = std::make_shared< FunctionType >( "xAgglomeration", agglomerationStorage_, level_, level_ );
      b_agglomeration_ = std::make_shared< FunctionType >( "bAgglomeration", agglomerationStorage_, level_, level_ );

      isStrategySet_ = true;
   }

   std::shared_ptr< PrimitiveStorage > originalStorage_;
   uint_t                              level_;
   bool                                solveOnEmptyProcesses_;
   bool                                isStrategySet_;
   bool                                isAgglomerationProcess_;

   std::shared_ptr< PrimitiveStorage > agglomerationStorage_;
   MigrationInfo                       migrationInfoToAgglomerationStorage_;
   MigrationInfo                       migrationInfoToOriginalStorage_;

   std::shared_ptr< OperatorType > A_agglomeration_;
   std::shared_ptr< FunctionType > b_agglomeration_;
   std::shared_ptr< FunctionType > x_agglomeration_;

   std::shared_ptr< Solver< OperatorType > > solver_;
};

} // namespace hyteg