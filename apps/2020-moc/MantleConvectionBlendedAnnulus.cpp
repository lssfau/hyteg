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

#include <cmath>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

/// In this benchmark we employ the Boussinesq-approximation and solve the constant-coefficient
/// Stokes eqn. coupled with velocity driven temperature transport.
///
/// The underlying domain is the blended annulus. At the inner domain boundary we apply heating
/// and cooling at the outer boundary. An initial temperature profile is given.
///
/// The convectivity of the non-dimensionalized setting is steered by the Rayleigh number.

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

std::shared_ptr< SetupPrimitiveStorage > buildSetupStorage( real_t rMin, real_t rMax, uint_t nTan, uint_t nRad )
{
   MeshInfo meshInfo     = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CRISS, nTan, nRad );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   AnnulusMap::setMap( *setupStorage );
   return setupStorage;
}

std::shared_ptr< Solver< P2P1ElementwiseBlendingStokesOperator > >
    buildStokesSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                       uint_t                                     minLevel,
                       uint_t                                     maxLevel,
                       uint_t                                     preSmooth,
                       uint_t                                     postSmooth,
                       real_t                                     uzawaOmega,
                       real_t                                     jacobiOmega )
{
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< P2P1ElementwiseBlendingStokesOperator, P1LumpedInvMassOperator > >(
           storage, minLevel, maxLevel );
   auto gaussSeidel = std::make_shared< WeightedJacobiSmoother< P2P1ElementwiseBlendingStokesOperator::VelocityOperator_T > >(
       storage, minLevel, maxLevel, jacobiOmega );
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< P2P1ElementwiseBlendingStokesOperator > >( storage,
                                                                                                                    gaussSeidel );
   auto smoother = std::make_shared< UzawaSmoother< P2P1ElementwiseBlendingStokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaOmega );
   auto coarseGridSolver = std::make_shared< MinResSolver< P2P1ElementwiseBlendingStokesOperator > >(
       storage, minLevel, minLevel, 1000, 1e-12, pressurePreconditioner );
   auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >();
   auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();

   auto gmgSolver = std::make_shared< GeometricMultigridSolver< P2P1ElementwiseBlendingStokesOperator > >( storage,
                                                                                                           smoother,
                                                                                                           coarseGridSolver,
                                                                                                           restrictionOperator,
                                                                                                           prolongationOperator,
                                                                                                           minLevel,
                                                                                                           maxLevel,
                                                                                                           preSmooth,
                                                                                                           postSmooth );
   return gmgSolver;
}

void runSimulation( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./MantleConvectionBlendedAnnulus.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   // Parameters

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t minLevel                = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel                = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t maxNumVCycles           = mainConf.getParameter< uint_t >( "maxNumVCycles" );
   const real_t stokesResidualTolerance = mainConf.getParameter< real_t >( "stokesResidualTolerance" );
   const uint_t stepsTotal              = mainConf.getParameter< uint_t >( "timesteps" );
   const real_t cflUpperBound           = mainConf.getParameter< real_t >( "cflUpperBound" );

   const uint_t preSmooth  = mainConf.getParameter< uint_t >( "preSmooth" );
   const uint_t postSmooth = mainConf.getParameter< uint_t >( "postSmooth" );

   const real_t rMin = mainConf.getParameter< real_t >( "rMin" );
   const real_t rMax = mainConf.getParameter< real_t >( "rMax" );
   const uint_t nTan = mainConf.getParameter< uint_t >( "nTan" );
   const uint_t nRad = mainConf.getParameter< uint_t >( "nRad" );

   // Should be > 0. The higher this number, the steeper the exp decay of the initial temp profile.
   const real_t initialTemperatureSteepness = mainConf.getParameter< real_t >( "initialTemperatureSteepness" );
   // Rayleigh number. With increasing Rayleigh number, the more does advection dominate diffusion.
   const real_t rayleighNumber = mainConf.getParameter< real_t >( "rayleighNumber" );

   const bool        writeVTK        = mainConf.getParameter< bool >( "writeVTK" );
   const bool        printTiming     = mainConf.getParameter< bool >( "printTiming" );
   const std::string outputDirectory = mainConf.getParameter< std::string >( "outputDirectory" );

   WALBERLA_LOG_INFO_ON_ROOT( " - level:        " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( " - time steps:   " << stepsTotal )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   // Domain

   auto setupStorage = buildSetupStorage( rMin, rMax, nTan, nRad );
   auto storage      = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   storage->getTimingTree()->start( "Total" );

   const real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( " - h (min, max): " << hMin << ", " << hMax )

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );

   // Functions and operators

   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > uLast( "uLast", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > outwardNormalField( "n", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel );

   P2Function< real_t > c( "c", storage, minLevel, maxLevel );

   P2P1ElementwiseBlendingStokesOperator L( storage, minLevel, maxLevel );
   P2ElementwiseBlendingMassOperator     M( storage, minLevel, maxLevel );

   const auto totalDoFsStokes = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "DoFs on max level: " << totalDoFsStokes );

   // Initial temperature

   std::function< real_t( const Point3D& ) > temperature = [&]( const Point3D& x ) {
      auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return std::exp( -initialTemperatureSteepness * ( ( radius - rMin ) / ( rMax - rMin ) ) );
   };

   c.interpolate( temperature, maxLevel );

   // Normal fields

   std::function< real_t( const hyteg::Point3D& ) > normalX = []( const hyteg::Point3D& x ) {
      return std::cos( std::atan2( x[1], x[0] ) );
   };

   std::function< real_t( const hyteg::Point3D& ) > normalY = []( const hyteg::Point3D& x ) {
      return std::sin( std::atan2( x[1], x[0] ) );
   };

   outwardNormalField.uvw().interpolate( { normalX, normalY }, maxLevel );

   // VTK

   if ( writeVTK )
   {
      writeDomainPartitioningVTK( storage, outputDirectory, "MantleConvectionBlendedAnnulusDomain" );
   }

   VTKOutput vtkOutput( outputDirectory, "MantleConvectionBlendedAnnulus", storage );

   vtkOutput.add( u );
   vtkOutput.add( f );
   vtkOutput.add( c );

   if ( writeVTK )
   {
      vtkOutput.write( maxLevel );
   }

   // Solver setup
#ifdef HYTEG_BUILD_WITH_PETSC
   WALBERLA_UNUSED( preSmooth );
   WALBERLA_UNUSED( postSmooth );
   auto stokesSolver = std::make_shared< PETScLUSolver< P2P1ElementwiseBlendingStokesOperator > >( storage, maxLevel );
#else
   auto stokesSolver = buildStokesSolver( storage, minLevel, maxLevel, preSmooth, postSmooth, 0.37, 0.66 );
#endif
   MMOCTransport< P2Function< real_t > > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   // Simulation loop

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       " %6s | %12s | %14s | %12s | %8s | %15s ", "step", "dt", "simulated time", "u max", "gmg iter", "Stokes residual" ) )

   real_t tTotal = 0;

   for ( uint_t i = 1; i <= stepsTotal; i++ )
   {
      uLast.assign( { 1.0 }, { u }, maxLevel, All );

      M.apply( c, f.uvw()[0], maxLevel, All );
      M.apply( c, f.uvw()[1], maxLevel, All );
      f.uvw().multElementwise( { f.uvw(), outwardNormalField.uvw() }, maxLevel );
      f.uvw().assign( { rayleighNumber }, { f.uvw() }, maxLevel, All );

      uint_t numVCycles      = 0;
      real_t currentResidual = std::numeric_limits< real_t >::max();
      while ( numVCycles < maxNumVCycles && currentResidual > stokesResidualTolerance )
      {
         stokesSolver->solve( L, u, f, maxLevel );
         L.apply( u, tmp, maxLevel, Inner | NeumannBoundary );
         r.assign( { 1.0, -1.0 }, { f, tmp }, maxLevel, Inner | NeumannBoundary );
         currentResidual = std::sqrt( r.dotGlobal( r, maxLevel, Inner | NeumannBoundary ) / real_c( totalDoFsStokes ) );
         numVCycles++;
         if ( !( numVCycles < maxNumVCycles && currentResidual > stokesResidualTolerance ) )
            break;
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " %6s | %12s | %14s | %12s | %8d | %15e ", "", "", "", "", numVCycles, currentResidual ) )
      }

      const auto maxVelocity = velocityMaxMagnitude( u.uvw(), tmp.uvw()[0], tmp.uvw()[1], maxLevel, All );
      const auto dt          = ( cflUpperBound / maxVelocity ) * hMin;

      transport.step( c, u.uvw(), uLast.uvw(), maxLevel, Inner, dt, 1, true );

      if ( writeVTK )
      {
         vtkOutput.write( maxLevel, i );
      }

      tTotal += dt;

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( " %6d | %12e | %14f | %12e | %8d | %15e ", i, dt, tTotal, maxVelocity, numVCycles, currentResidual ) )
   }

   storage->getTimingTree()->stop( "Total" );
   if ( printTiming )
   {
      WALBERLA_LOG_INFO_ON_ROOT( storage->getTimingTree()->getCopyWithRemainder() );
   }
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::runSimulation( argc, argv );
   return EXIT_SUCCESS;
}
