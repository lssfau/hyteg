/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include <hyteg/solvers/GaussSeidelSmoother.hpp>
#include <hyteg/solvers/controlflow/SolverLoop.hpp>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/Git.hpp"
#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1toP1InjectionRestriction.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

using walberla::real_c;
using walberla::real_t;

namespace hyteg {

void simulate( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   printGitInfo();
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./StokesSphereTransport.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   /////////////// Parameters ///////////////
   const walberla::Config::BlockHandle mainConf    = cfg->getBlock( "Parameters" );
   const walberla::Config::BlockHandle layersParam = cfg->getBlock( "Layers" );

   const uint_t          ntan = mainConf.getParameter< uint_t >( "ntan" );
   std::vector< double > layers;
   for ( const auto& it : layersParam )
   {
      layers.push_back( layersParam.getParameter< double >( it.first ) );
   }

   const double      rmin                 = layers.front();
   const double      rmax                 = layers.back();
   const uint_t      minLevel             = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel             = mainConf.getParameter< uint_t >( "maxLevel" );
   const real_t      stokesResidual       = mainConf.getParameter< real_t >( "stokesResidual" );
   const uint_t      stokesMaxNumVCycles  = mainConf.getParameter< uint_t >( "stokesMaxNumVCycles" );
   const uint_t      stokesSolveInterval  = mainConf.getParameter< uint_t >( "stokesSolveInterval" );
   const uint_t      numDiffusionVCycles  = mainConf.getParameter< uint_t >( "numDiffusionVCycles" );
   const uint_t      timeSteps            = mainConf.getParameter< uint_t >( "timeSteps" );
   const real_t      diffusivity          = mainConf.getParameter< real_t >( "diffusivity" );
   const real_t      dt                   = mainConf.getParameter< real_t >( "dt" );
   const real_t      rhsScaleFactor       = mainConf.getParameter< real_t >( "rhsScaleFactor" );
   const bool        writeVTK             = mainConf.getParameter< bool >( "vtkOutput" );
   const bool        writeDomainVTK       = mainConf.getParameter< bool >( "writeDomainVTK" );
   const bool        exitAfterWriteDomain = mainConf.getParameter< bool >( "exitAfterWriteDomain" );
   const uint_t      VTKOutputFrequency   = mainConf.getParameter< uint_t >( "vtkFrequency" );
   const bool        printTiming          = mainConf.getParameter< bool >( "printTiming" );
   const std::string timingFile           = mainConf.getParameter< std::string >( "timingFile" );
   const std::string vtkBaseFile          = mainConf.getParameter< std::string >( "vtkBaseFile" );
   const std::string vtkDirectory         = mainConf.getParameter< std::string >( "vtkDirectory" );
   const uint_t      vtkOutputLevel       = mainConf.getParameter< uint_t >( "vtkOutputLevel" );

   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" )
   WALBERLA_LOG_INFO_ON_ROOT( " - domain:" )
   WALBERLA_LOG_INFO_ON_ROOT( "   + rmin: " << rmin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + rmax: " << rmax )
   WALBERLA_LOG_INFO_ON_ROOT( "   + layers (min == 2): " << layers.size() )
   WALBERLA_LOG_INFO_ON_ROOT( "   + ntan: " << ntan )
   WALBERLA_LOG_INFO_ON_ROOT( " - simulation:" )
   WALBERLA_LOG_INFO_ON_ROOT( "   + diffusivity: " << diffusivity )
   WALBERLA_LOG_INFO_ON_ROOT( "   + rhs scale factor: " << rhsScaleFactor )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dt: " << dt )
   WALBERLA_LOG_INFO_ON_ROOT( " - solvers:" )
   WALBERLA_LOG_INFO_ON_ROOT( "   + min level: " << minLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + max level: " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Stokes target residual: " << stokesResidual )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Stokes max v-cycles: " << stokesMaxNumVCycles )
   WALBERLA_LOG_INFO_ON_ROOT( "   + diffusion v-cycles: " << numDiffusionVCycles )
   WALBERLA_LOG_INFO_ON_ROOT( " - other:" )
   WALBERLA_LOG_INFO_ON_ROOT( "   + write VTK: " << writeVTK )
   WALBERLA_LOG_INFO_ON_ROOT( "   + write domain: " << writeDomainVTK )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK directory: " << vtkDirectory )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK base name: " << vtkBaseFile )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK output level: " << vtkOutputLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + VTK interval: " << VTKOutputFrequency )
   WALBERLA_LOG_INFO_ON_ROOT( "   + print timing: " << printTiming )
   WALBERLA_LOG_INFO_ON_ROOT( "   + exit after domain output: " << exitAfterWriteDomain )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   /////////////////// Mesh / Domain ///////////////////////

   MeshInfo meshInfo     = MeshInfo::meshSphericalShell( ntan, layers );
   auto     setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< PrimitiveStorage >( *setupStorage, timingTree, 1 );

   storage->getTimingTree()->start( "Total" );

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );

   if ( writeDomainVTK )
   {
      writeDomainPartitioningVTK( storage, vtkDirectory, vtkBaseFile + "_domain" );
   }

   const auto hmin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const auto hmax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "hmin: " << hmin )
   WALBERLA_LOG_INFO_ON_ROOT( "hmax: " << hmax );
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   if ( exitAfterWriteDomain )
      return;

   P2P1TaylorHoodFunction< real_t > r( "r", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > uLastTimeStep( "uLast", storage, minLevel, maxLevel );
   P2Function< real_t >             temp( "temperature", storage, minLevel, maxLevel );
   P2Function< real_t >             tempOld( "temperature_old", storage, minLevel, maxLevel );
   P2Function< real_t >             tempTmp( "temperature_tmp", storage, minLevel, maxLevel );
   P2Function< real_t >             tempR( "temperature_residual", storage, minLevel, maxLevel );
   P2Function< real_t >             normalX( "normalX", storage, minLevel, maxLevel );
   P2Function< real_t >             normalY( "normalY", storage, minLevel, maxLevel );
   P2Function< real_t >             normalZ( "normalZ", storage, minLevel, maxLevel );
   P2Function< real_t >             tmp( "tmp", storage, minLevel, maxLevel );
   P2Function< real_t >             tmp2( "tmp2", storage, minLevel, maxLevel );

   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      uint_t tmpDofStokes = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "Stokes DoFs on level " << lvl << " : " << tmpDofStokes );
      uint_t tmpDofTemperature = numberOfGlobalDoFs< P2FunctionTag >( *storage, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "Temperature DoFs on level " << lvl << " : " << tmpDofTemperature );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   VTKOutput vtkOutput( vtkDirectory, vtkBaseFile, storage, VTKOutputFrequency );
   if ( writeVTK )
   {
      // vtkOutput.add( u.u.getVertexDoFFunction() );
      // vtkOutput.add( u.v.getVertexDoFFunction() );
      // vtkOutput.add( u.w.getVertexDoFFunction() );
      vtkOutput.add( temp.getVertexDoFFunction() );
   }
   P1toP1InjectionRestriction injectionRestriction;

   P2P1TaylorHoodStokesOperator L( storage, minLevel, maxLevel );
   P2ConstantLaplaceOperator    laplace( storage, minLevel, maxLevel );
   P2ConstantMassOperator       M( storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > temperature = [rmin, rmax]( const Point3D& x ) {
      return std::pow( ( rmax - x.norm() ) / ( rmax - rmin ), 3.0 );
   };

   temp.interpolate( temperature, maxLevel );

   std::function< real_t( const Point3D& ) > zero = []( const Point3D& ) { return 0.0; };
   std::function< real_t( const Point3D& ) > ones = []( const Point3D& ) { return 1.0; };

   std::function< real_t( const Point3D& ) > nX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > nY = []( const Point3D& x ) { return x[1] / x.norm(); };
   std::function< real_t( const Point3D& ) > nZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   normalX.interpolate( nX, maxLevel );
   normalY.interpolate( nY, maxLevel );
   normalZ.interpolate( nZ, maxLevel );

   auto coarseGridSolver = std::make_shared< PETScBlockPreconditionedStokesSolver< P2P1TaylorHoodStokesOperator > >(
       storage, minLevel, 1e-12, 2000, 1 );
   auto stokesRestriction  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );
   auto stokesProlongation = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto gaussSeidel        = std::make_shared< hyteg::GaussSeidelSmoother< P2P1TaylorHoodStokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner =
       std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< P2P1TaylorHoodStokesOperator > >( storage,
                                                                                                                  gaussSeidel );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< P2P1TaylorHoodStokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, 0.3 );
   const auto uzawaRelaxationParameter = estimateUzawaRelaxationParameter( storage, uzawaVelocityPreconditioner, 2, 20, 2 );
   uzawaSmoother->setRelaxationParameter( uzawaRelaxationParameter );
   WALBERLA_LOG_INFO_ON_ROOT( "Estimated relaxation parameter for Uzawa: " << uzawaRelaxationParameter );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   auto gmgSolver = std::make_shared< GeometricMultigridSolver< P2P1TaylorHoodStokesOperator > >(
       storage, uzawaSmoother, coarseGridSolver, stokesRestriction, stokesProlongation, minLevel, maxLevel, 3, 3, 2 );

   MMOCTransport< P2Function< real_t > > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );

   P2ConstantUnsteadyDiffusionOperator diffusionOperator(
       storage, minLevel, maxLevel, dt, diffusivity, DiffusionTimeIntegrator::ImplicitEuler );
   auto diffusionCoarseGridSolver =
       std::make_shared< CGSolver< P2ConstantUnsteadyDiffusionOperator > >( storage, minLevel, maxLevel, 2000, 1e-12 );
   auto diffusionRestriction  = std::make_shared< P2toP2QuadraticRestriction >();
   auto diffusionProlongation = std::make_shared< P2toP2QuadraticProlongation >();
   auto gsSmoother            = std::make_shared< GaussSeidelSmoother< P2ConstantUnsteadyDiffusionOperator > >();
   auto diffusionGMGSolver    = std::make_shared< GeometricMultigridSolver< P2ConstantUnsteadyDiffusionOperator > >(
       storage, gsSmoother, diffusionCoarseGridSolver, diffusionRestriction, diffusionProlongation, minLevel, maxLevel, 2, 1, 0 );
   auto diffusionSolver =
       std::make_shared< SolverLoop< P2ConstantUnsteadyDiffusionOperator > >( diffusionGMGSolver, numDiffusionVCycles );

   UnsteadyDiffusion< P2Function< real_t >,
                      P2ConstantUnsteadyDiffusionOperator,
                      P2ConstantLaplaceOperator,
                      P2ConstantMassOperator >
       diffusion( storage, minLevel, maxLevel, diffusionSolver );

   printFunctionAllocationInfo( *storage, 1 );

   walberla::WcTimer timer;

   auto calculateResidualStokes = [&]() {
      L.apply( u, r, maxLevel, Inner | NeumannBoundary );
      r.assign( { 1.0, -1.0 }, { f, r }, maxLevel, Inner | NeumannBoundary );

      real_t res = sqrt( r.dotGlobal( r, maxLevel, All ) ) /
                   real_c( numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel ) );
      return res;
   };

   auto calculateResidualDiffusion = [&]() {
      diffusionOperator.apply( temp, tempR, maxLevel, Inner | NeumannBoundary );
      M.apply( tempOld, tempTmp, maxLevel, Inner | NeumannBoundary );
      tempR.assign( { 1.0, -1.0 }, { tempTmp, tempR }, maxLevel, Inner | NeumannBoundary );

      real_t res =
          sqrt( tempR.dotGlobal( tempR, maxLevel, All ) ) / real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel ) );
      return res;
   };

   auto writeVTKCallback = [&]( uint_t timestep ) {
      storage->getTimingTree()->start( "VTK" );

      WALBERLA_LOG_INFO_ON_ROOT( "VTK output ..." )
      timer.start();

      for ( uint_t sourceLevel = maxLevel; sourceLevel > vtkOutputLevel; sourceLevel-- )
      {
         injectionRestriction.restrict( temp.getVertexDoFFunction(), sourceLevel, All );
      }

      vtkOutput.write( vtkOutputLevel, timestep );

      timer.end();
      WALBERLA_LOG_INFO_ON_ROOT( "" )
      WALBERLA_LOG_INFO_ON_ROOT( "... " << timer.last() << " seconds" )
      WALBERLA_LOG_INFO_ON_ROOT( "" )

      storage->getTimingTree()->stop( "VTK" );
   };

   auto maxMagnitudeVelocity = [&]() {
      tmp2.interpolate( 0, maxLevel, All );

      tmp.assign( { 1.0 }, { u.uvw()[0] }, maxLevel, All );
      tmp.multElementwise( { tmp, tmp }, maxLevel, All );
      tmp2.assign( { 1.0, 1.0 }, { tmp2, tmp }, maxLevel, All );

      tmp.assign( { 1.0 }, { u.uvw()[1] }, maxLevel, All );
      tmp.multElementwise( { tmp, tmp }, maxLevel, All );
      tmp2.assign( { 1.0, 1.0 }, { tmp2, tmp }, maxLevel, All );

      tmp.assign( { 1.0 }, { u.uvw()[2] }, maxLevel, All );
      tmp.multElementwise( { tmp, tmp }, maxLevel, All );
      tmp2.assign( { 1.0, 1.0 }, { tmp2, tmp }, maxLevel, All );

      return std::sqrt( tmp2.getMaxMagnitude( maxLevel, All ) );
   };

   if ( writeVTK )
   {
      writeVTKCallback( 0 );
   }

   storage->getTimingTree()->start( "Simulation" );

   for ( uint_t step = 0; step < timeSteps; ++step )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "##### Time step " << step << " #####" )
      WALBERLA_LOG_INFO_ON_ROOT( "" )

      uLastTimeStep.assign( { 1.0 }, { u }, maxLevel, All );

      // Updating right-hand side (Boussinesq approximation)

      M.apply( temp, f.uvw()[0], maxLevel, All );
      M.apply( temp, f.uvw()[1], maxLevel, All );
      M.apply( temp, f.uvw()[2], maxLevel, All );

      f.uvw()[0].multElementwise( { f.uvw()[0], normalX }, maxLevel, All );
      f.uvw()[1].multElementwise( { f.uvw()[1], normalY }, maxLevel, All );
      f.uvw()[2].multElementwise( { f.uvw()[2], normalZ }, maxLevel, All );

      f.uvw().assign( { rhsScaleFactor }, { f.uvw() }, maxLevel, All );

      // Stokes solver

      // Check residual, if less than tolerance, quit v-cycle loop early

      if ( step % stokesSolveInterval == 0 )
      {
         storage->getTimingTree()->start( "Stokes" );
         timer.start();

         real_t currentResidual = calculateResidualStokes();
         real_t lastResidual    = currentResidual;
         WALBERLA_LOG_INFO_ON_ROOT( " iteration | residual (l2) | convergence rate " );
         WALBERLA_LOG_INFO_ON_ROOT( "-----------+---------------+------------------" );
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "   initial | %13.3e |                - ", currentResidual ) )

         for ( uint_t i = 0; i < stokesMaxNumVCycles; i++ )
         {
            if ( currentResidual < stokesResidual )
               break;

            gmgSolver->solve( L, u, f, maxLevel );

            currentResidual      = calculateResidualStokes();
            auto convergenceRate = currentResidual / lastResidual;
            lastResidual         = currentResidual;
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %9d | %13.3e | %16.3e ", i, currentResidual, convergenceRate ) )
         }

         const auto maxVelocity = maxMagnitudeVelocity();
         const auto cfl         = ( maxVelocity * dt ) / hmin;
         WALBERLA_LOG_INFO_ON_ROOT( "" )
         WALBERLA_LOG_INFO_ON_ROOT( "max velocity magnitude: " << maxVelocity );
         WALBERLA_LOG_INFO_ON_ROOT( "CFL ( u_max * dt / h_min ): " << cfl );

         timer.end();
         WALBERLA_LOG_INFO_ON_ROOT( "" )
         WALBERLA_LOG_INFO_ON_ROOT( "... " << timer.last() << " seconds" )
         WALBERLA_LOG_INFO_ON_ROOT( "" )
         storage->getTimingTree()->stop( "Stokes" );
      }

      // Advection-diffusion

      storage->getTimingTree()->start( "Advection" );

      WALBERLA_LOG_INFO_ON_ROOT( "Advection time step ..." )

      timer.start();

      transport.step( temp, u.uvw(), uLastTimeStep.uvw(), maxLevel, All, dt, 1, true );

      timer.end();
      WALBERLA_LOG_INFO_ON_ROOT( "" )
      WALBERLA_LOG_INFO_ON_ROOT( "... " << timer.last() << " seconds" )
      WALBERLA_LOG_INFO_ON_ROOT( "" )

      storage->getTimingTree()->stop( "Advection" );

      if ( diffusivity > 0 )
      {
         storage->getTimingTree()->start( "Diffusion" );

         WALBERLA_LOG_INFO_ON_ROOT( "Diffusion time step ..." )

         timer.start();

         tempOld.assign( { 1.0 }, { temp }, maxLevel, All );
         diffusion.step( diffusionOperator, laplace, M, temp, tempOld, maxLevel, Inner | NeumannBoundary );

         const auto residualDiffusion = calculateResidualDiffusion();
         WALBERLA_LOG_INFO_ON_ROOT( "" )
         WALBERLA_LOG_INFO_ON_ROOT( "l2 residual diffusion: " << residualDiffusion );

         timer.end();
         WALBERLA_LOG_INFO_ON_ROOT( "" )
         WALBERLA_LOG_INFO_ON_ROOT( "... " << timer.last() << " seconds" )
         WALBERLA_LOG_INFO_ON_ROOT( "" )

         storage->getTimingTree()->stop( "Diffusion" );
      }

      if ( writeVTK )
      {
         writeVTKCallback( step + 1 );
      }
   }

   storage->getTimingTree()->stop( "Simulation" );

   storage->getTimingTree()->stop( "Total" );

   if ( printTiming )
   {
      printTimingTree( *storage->getTimingTree() );
      writeTimingTreeJSON( *storage->getTimingTree(), vtkDirectory + "/" + timingFile );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   hyteg::simulate( argc, argv );
}
