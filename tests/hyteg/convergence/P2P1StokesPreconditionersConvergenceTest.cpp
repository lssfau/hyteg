/*
* Copyright (c) 2017-2024 Andreas Burkhart, Nils Kohl.
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

#include "core/math/Random.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/SubstitutePreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"
#include "hyteg_operators_composites/viscousblock/P2ViscousBlockLaplaceOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

// Tests block-diagonal Stokes preconditioners with the generated operators.
//
// The test case checks the convergence of the solver for the "full" Stokes operator and a mildly varying viscosity.
// It's certainly not the optimal choice of components and parameters, just there to test the functionality.
//
// It employs the "substitute" operators that facilitate plugging together approximate preconditioners.
// They enable, for instance, the approximation of the inverse full viscous block with a few multigrid iterations on the
// block-Laplace operator. Overall, the solver looks like the following:
//
// Outer solver:   FGMRES
// Preconditioner: inexact Uzawa iteration
//                 - A-block approximation in the inexact Uzawa: multigrid-preconditioned CG
//                   (the multigrid only acts on the block-Laplace operator)
//                 - Schur complement approximation in the inexact Uzawa: CG on inversely-viscosity weighted pressure-mass matrix
//

template < typename StokesFullOperator, typename BlockLaplaceOperator, typename SchurOperator >
void testFullViscousStokes( uint_t dim, uint_t maxLevel, uint_t fgmresIterations, real_t unscaledResidualEpsilon )
{
   const bool   writeVTK                          = false;
   const uint_t minLevel                          = 2;
   const uint_t chebyshevSpectralRadiusIterations = 20;

   walberla::math::seedRandomGenerator( 42 );

   auto rand = []( const Point3D& ) { return walberla::math::realRandom(); };
   auto one  = []( const Point3D& ) { return real_c( 1.0 ); };

   const auto blending =
       std::is_same_v< StokesFullOperator, operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator > && dim == 3;

   std::shared_ptr< PrimitiveStorage > storage;

   if ( blending )
   {
      WALBERLA_CHECK_EQUAL( dim, 3, "Blending test only implemented for 3D (IcoShell)." );
      MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 3, 2, 0.5, 1.0 );
      SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      IcosahedralShellMap::setMap( setupStorage );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }
   else
   {
      WALBERLA_CHECK_EQUAL( dim, 2, "Non-blending test only implemented for 2D." );
      MeshInfo              meshInfo = MeshInfo::meshUnitSquare( 0 );
      SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      storage = std::make_shared< PrimitiveStorage >( setupStorage );
   }

   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );

   P2P1TaylorHoodFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   P2Function< real_t > mu( "mu", storage, minLevel, maxLevel );
   P1Function< real_t > muInv( "muInv", storage, minLevel, maxLevel );

   // Some coefficient and its inverse.
   auto visc    = [&]( const Point3D& x ) { return 2 + std::sin( x[0] ) * std::cos( x[1] ); };
   auto viscInv = [&]( const Point3D& x ) { return 1.0 / visc( x ); };
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      mu.interpolate( visc, level, All );
      muInv.interpolate( viscInv, level, All );
   }

   // Evaluating linear form for right-hand side.
   L2Space< 7, P2Function< real_t >, real_t > l2space( storage, maxLevel );
   l2space.dot( one, f.uvw()[0] );
   l2space.dot( one, f.uvw()[1] );
   if ( storage->hasGlobalCells() )
   {
      l2space.dot( one, f.uvw()[2] );
   }

   // Operators
   using SubstAType = BlockLaplaceOperator;
   StokesFullOperator L( storage, minLevel, maxLevel, mu );
   SchurOperator      S( storage, minLevel, maxLevel, muInv );

   L.getA().computeInverseDiagonalOperatorValues();

   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "------- Setup solvers & operators: Starting -------" );
   WALBERLA_LOG_INFO_ON_ROOT( "---------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "Setup Stokes prolongation and restriction" );

   WALBERLA_LOG_INFO_ON_ROOT( "Setup approximate A-block solver" );

   auto APrecOperator = std::make_shared< SubstAType >( storage, minLevel, maxLevel );
   APrecOperator->computeInverseDiagonalOperatorValues();

   auto ABlockProlongationOperator = std::make_shared< P2toP2QuadraticVectorProlongation >();
   auto ABlockRestrictionOperator  = std::make_shared< P2toP2QuadraticVectorRestriction >();
   auto ABlockCoarseGridSolver     = std::make_shared< CGSolver< SubstAType > >( storage, minLevel, maxLevel, 10, 1e-8 );
   auto ABlockSmoother             = std::make_shared< ChebyshevSmoother< SubstAType > >( storage, minLevel, maxLevel );

   // Avoiding startpoint of the power iteration to be within the operator kernel
   tmp.uvw().interpolate( rand, maxLevel, All );
   auto spectralRadiusA =
       chebyshev::estimateRadius( *APrecOperator, maxLevel, chebyshevSpectralRadiusIterations, storage, tmp.uvw(), r.uvw() );

   ABlockSmoother->setupCoefficients( 3, spectralRadiusA );
   auto ABlockMultigridSolver = std::make_shared< GeometricMultigridSolver< SubstAType > >( storage,
                                                                                            ABlockSmoother,
                                                                                            ABlockCoarseGridSolver,
                                                                                            ABlockRestrictionOperator,
                                                                                            ABlockProlongationOperator,
                                                                                            minLevel,
                                                                                            maxLevel,
                                                                                            1,
                                                                                            1,
                                                                                            0,
                                                                                            CycleType::VCYCLE );

   auto SubstABlockSolver =
       std::make_shared< SubstitutePreconditioner< typename StokesFullOperator::ViscousOperator_T, SubstAType > >(
           ABlockMultigridSolver, APrecOperator );

   auto ABlockSolver = std::make_shared< CGSolver< typename StokesFullOperator::ViscousOperator_T > >(
       storage, minLevel, maxLevel, 3, 1e-2, SubstABlockSolver );

   WALBERLA_LOG_INFO_ON_ROOT( "Setup approx. Schur complement solver" );

   auto SchurSolver = std::make_shared< CGSolver< SchurOperator > >( storage, minLevel, maxLevel, 50, 1e-8 );

   WALBERLA_LOG_INFO_ON_ROOT( "Setup Stokes preconditioner" );

   uint_t chebyshevIterationsPerUzawa = 1;

   auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >();

   auto uzawaSmoother = std::make_shared<
       InexactUzawaPreconditioner< StokesFullOperator, typename StokesFullOperator::ViscousOperator_T, SchurOperator > >(
       storage, minLevel, maxLevel, S, ABlockSolver, SchurSolver, 1.0, 1.0, chebyshevIterationsPerUzawa );

   auto coarseGridSolver = solvertemplates::stokesMinResSolver< StokesFullOperator >( storage, minLevel, 1e-8, 10 );

   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesFullOperator > >( storage,
                                                                                              uzawaSmoother,
                                                                                              coarseGridSolver,
                                                                                              restrictionOperator,
                                                                                              prolongationOperator,
                                                                                              minLevel,
                                                                                              maxLevel,
                                                                                              2,
                                                                                              2,
                                                                                              2,
                                                                                              CycleType::VCYCLE );

   WALBERLA_LOG_INFO_ON_ROOT( "Setup Stokes solver" );

   auto solver = std::make_shared< FGMRESSolver< StokesFullOperator > >(
       storage, minLevel, maxLevel, fgmresIterations, 50, 1e-8, 1e-8, 0, multigridSolver );
   solver->setPrintInfo( true );

   WALBERLA_LOG_INFO_ON_ROOT( "Starting FGMRES" )

   // Initial residual.
   L.apply( u, r, maxLevel, Inner | NeumannBoundary );
   r.assign( { 1.0, -1.0 }, { r, f }, maxLevel );
   real_t initialUnscaledResiduum = std::sqrt( r.dotGlobal( r, maxLevel, Inner ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Initial residual: " << initialUnscaledResiduum )

   solver->solve( L, u, f, maxLevel );

   L.apply( u, r, maxLevel, Inner | NeumannBoundary );
   r.assign( { 1.0, -1.0 }, { r, f }, maxLevel );

   real_t unscaledFinalResiduum = std::sqrt( r.dotGlobal( r, maxLevel, Inner ) );

   WALBERLA_LOG_INFO_ON_ROOT( "Final residual: " << unscaledFinalResiduum )
   WALBERLA_CHECK_LESS( unscaledFinalResiduum, unscaledResidualEpsilon );

   // printTimingTree( *storage->getTimingTree() );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output", "StokesConvTest", storage );
      vtk.add( mu );
      vtk.add( u );
      vtk.write( maxLevel );
   }
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   bool longrun = false;
   for ( int i = 0; i < argc; i++ )
   {
      auto arg = std::string( argv[i] );
      if ( arg == "--longrun" )
      {
         longrun = true;
      }
   }

   {
      using StokesFullOperator   = operatorgeneration::P2P1StokesFullOperator;
      using BlockLaplaceOperator = operatorgeneration::P2ViscousBlockLaplaceOperator;
      using SchurOperator        = operatorgeneration::P1ElementwiseKMass;

      testFullViscousStokes< StokesFullOperator, BlockLaplaceOperator, SchurOperator >( 2, 4, 10, 1e-8 );
   }

   if ( longrun )
   {
      using StokesFullOperator   = operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator;
      using BlockLaplaceOperator = operatorgeneration::P2ViscousBlockLaplaceIcosahedralShellMapOperator;
      using SchurOperator        = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

      testFullViscousStokes< StokesFullOperator, BlockLaplaceOperator, SchurOperator >( 3, 3, 5, 1.5e-4 );
   }

   return EXIT_SUCCESS;
}
