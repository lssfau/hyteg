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
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/SymmetricGaussSeidelSmoother.hpp"
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

template < typename config >
static void solvePoisson( uint_t minLevel, uint_t maxLevel, uint_t cycles )
{
   const uint_t numEdgesPerSide = 1;
   auto         meshInfo =
       MeshInfo::meshCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), numEdgesPerSide, numEdgesPerSide, numEdgesPerSide );

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
   auto smoother             = std::make_shared< typename config::smooth >();
   auto restrictionOperator  = std::make_shared< typename config::Restriction >();
   auto prolongationOperator = std::make_shared< typename config::Prolongation >();
   auto gmgSolver            = std::make_shared< GeometricMultigridSolver< typename config::laplaceOP > >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   std::map< std::string, std::function< void( uint_t ) > > differentSolvers;
   std::map< std::string, double >                          timings;

   timer.reset();

   for ( int i = 0; i < cycles; ++i )
   {
      gmgSolver->solve( A, u, f, maxLevel );
   }
   timer.end();

   WALBERLA_LOG_INFO( walberla::format( "Max Value: %9.3e | Time: %9.3e", u.getMaxValue( maxLevel, Inner ), timer.last() ) );
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager;
   const uint_t        minLevel = 2;
   const uint_t        maxLevel = 4;
   using P1op                   = hyteg::P1ConstantLaplaceOperator;
   using P2op                   = hyteg::P2ConstantLaplaceOperator;

   hyteg::solvePoisson< hyteg::P1conf< hyteg::CGSolver< P1op >, hyteg::GaussSeidelSmoother< P1op > > >( minLevel, maxLevel, 1 );
   hyteg::solvePoisson< hyteg::P1conf< hyteg::MinResSolver< P1op >, hyteg::SymmetricGaussSeidelSmoother< P1op > > >(
       minLevel, maxLevel, 2 );

   hyteg::solvePoisson< hyteg::P2conf< hyteg::CGSolver< P2op >, hyteg::GaussSeidelSmoother< P2op > > >( minLevel, maxLevel, 1 );
   hyteg::solvePoisson< hyteg::P2conf< hyteg::MinResSolver< P2op >, hyteg::SymmetricGaussSeidelSmoother< P2op > > >(
       minLevel, maxLevel, 1 );
   return 0;
}