/*
* Copyright (c) 2017-2023 Dominik Thoennes.
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

#include "core/Environment.h"
#include "core/extern/json.hpp"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/SORSmoother.hpp"
#include "hyteg/solvers/SymmetricSORSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"

using walberla::real_c;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

template < class coarseSolver, class smoother >
struct P1conf
{
   using function     = P1Function< real_t >;
   using laplaceOP    = P1ConstantLaplaceOperator;
   using Prolongation = P1toP1LinearProlongation<>;
   using Restriction  = P1toP1LinearRestriction< real_t >;
   using cs           = coarseSolver;
   using smooth       = smoother;
};

template < class coarseSolver, class smoother >
struct P2conf
{
   using function     = P2Function< real_t >;
   using laplaceOP    = P2ConstantLaplaceOperator;
   using Prolongation = P2toP2QuadraticProlongation;
   using Restriction  = P2toP2QuadraticRestriction;
   using dt           = real_t;
   using cs           = coarseSolver;
   using smooth       = smoother;
};

///
/// \tparam config can be "P1conf" or "P2conf"
/// \param minLevel the lowest level in the Multigrid hierarchy; should be >= 2
/// \param maxLevel the highest level in the Multigird hierarchy; must be >= minLevel
/// \param cycles numer of Multigrid cycles; must be >= 1
/// \param smoothParameter Parameter for the smoother; must be 2 > x > 1 for SOR and 1 > x > 0 for weighted Jacobi
template < typename config >
static void solvePoisson( uint_t minLevel, uint_t maxLevel, uint_t cycles, real_t smoothParameter )
{
   const uint_t numEdgesPerSide   = 1;
   const uint_t coarseRefinements = 0;

   auto meshInfo =
       MeshInfo::meshCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), numEdgesPerSide, numEdgesPerSide, numEdgesPerSide );
   meshInfo = MeshInfo::refinedCoarseMesh( meshInfo, coarseRefinements );

   std::function< real_t( const hyteg::Point3D& ) > random = []( const hyteg::Point3D& ) { return walberla::math::realRandom(); };

   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   typename config::function u( "u", storage, minLevel, maxLevel );
   typename config::function f( "f", storage, minLevel, maxLevel );

   walberla::timing::Timer< walberla::timing::WcPolicy > timer;

   typename config::laplaceOP A( storage, minLevel, maxLevel );

   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      u.interpolate( 0, l, DirichletBoundary );
      u.interpolate( random, l, Inner );
      f.setToZero( l );
   }

   auto coarseGridSolver     = std::make_shared< typename config::cs >( storage, minLevel, maxLevel );
   auto smoother             = std::make_shared< typename config::smooth >( smoothParameter );
   auto restrictionOperator  = std::make_shared< typename config::Restriction >();
   auto prolongationOperator = std::make_shared< typename config::Prolongation >();
   auto gmgSolver            = std::make_shared< GeometricMultigridSolver< typename config::laplaceOP > >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   std::map< std::string, std::function< void( uint_t ) > > differentSolvers;
   std::map< std::string, double >                          timings;

   timer.reset();

   for ( uint_t i = 0; i < cycles; ++i )
   {
      gmgSolver->solve( A, u, f, maxLevel );
   }
   timer.end();

   WALBERLA_LOG_INFO(
       walberla::format( "Max Magnitude: %9.3e | Time: %9.3e", u.getMaxMagnitude( maxLevel, Inner ), timer.last() ) );
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   const uint_t        minLevel = 2;
   const uint_t        maxLevel = 4;
   using P1op                   = hyteg::P1ConstantLaplaceOperator;
   using P2op                   = hyteg::P2ConstantLaplaceOperator;

   // P1conf and P2conf require two template parameter
   // - coarseSolver is the solver that is used on the coarsest grid in the Multigrid hierarchy; this can be:
   //   CGSolver, MinResSolver or GMRESSolver
   // - smoother is the smoother that is used on each level in the smoother Multigrid hierarchy; this can be:
   //   SORSmoother or SymmetricSORSmoother
   //   the symmetric variant alternates the direction of the smoothing for each iteration

   hyteg::solvePoisson< hyteg::P1conf< hyteg::CGSolver< P1op >, hyteg::SORSmoother< P1op > > >( minLevel, maxLevel, 1, 1.0 );
   hyteg::solvePoisson< hyteg::P1conf< hyteg::MinResSolver< P1op >, hyteg::SymmetricSORSmoother< P1op > > >(
       minLevel, maxLevel, 2, 0.3 );

   hyteg::solvePoisson< hyteg::P2conf< hyteg::CGSolver< P2op >, hyteg::SORSmoother< P2op > > >( minLevel, maxLevel, 1, 0.5 );
   hyteg::solvePoisson< hyteg::P2conf< hyteg::GMRESSolver< P2op >, hyteg::SORSmoother< P2op > > >( minLevel, maxLevel, 1, 0.3 );
   return 0;
}