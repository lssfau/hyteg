/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

/**
 * \page 06_FullAppPlumeInCube A full app creating a plume in a cube
 *
 * \dontinclude tutorials/06_FullAppPlumeInCube/Plume.cpp
 *
 * \brief In this tutorial we will set up a complete app that will solve a coupled system
 * of the Stokes equations and the convection-diffusion equation.
 *
 * \section FullAppPlumeInCube-param Parameter file reader
 *
 * As in all applications, we first setup the MPI environment.
 * We also load a parameter file to change parameters without having to rebuild the application.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Setup environment
 *
 * \section FullAppPlumeInCube-domain Domain
 *
 * Next, we define our domain. In this tutorial, we will use one of the internal
 * mesh generators to create a cuboid mesh.
 *
 * From this we define our SetupPrimitiveStorage and set the boundary flags.
 * A boundary flag is attached to each primitive - we use a convenience function
 * to set all flags of all primitives on the boundary to 1 and all others to 0.
 * Later, this helps us to define boundary conditions.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Domain
 *
 * \note
 * Each individual function variable may treat the boundary flags at the primitives differently.
 * However, mostly (also in this app) the defaults are convenient:
 * - 0 -> inner domain
 * - 1 -> Dirichlet boundary
 * - 2 -> Neumann boundary
 *
 * \note
 * For the Stokes functions the default is slightly different: the pressure variables are treating
 * all flags as inner domain. This makes sense since we do not really have boundary conditions for
 * the pressure in theory.
 *
 * \section FullAppPlumeInCube-discretization Discretization, function spaces and operators
 *
 * Now we define the function spaces and respective variables as well as our discrete operators.
 * Alternately, we solve the Stokes equation to obtain a velocity and pressure field from an
 * external force. Then we use the resulting velocity field to update the transport direction
 * of the temperature field. Employing a Boussinesq approximation we model buoyancy directly
 * through the temperature in z-direction. Therefore we plug in the temperature as the right-hand
 * side of the momentum equation.
 *
 * To advance in time we employ an explicit, algebraic upwind scheme. The Stokes equation is solved every few time
 * steps using a monolithic multigrid method with inexact-Uzawa type relaxation. All quantities are discretized
 * with linear finite elements. Therefore we employ a PSPG stabilization for the pressure.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Discretization
 *
 * \section FullAppPlumeInCube-bc Initial and boundary conditions
 *
 * Now we interpolate the initial conditions. Boundary conditions for the velocity are no-slip
 * (== Dirichlet, u = 0) everywhere. Since this is the default, we do not need do anything here.
 *
 * For the temperature, we choose to interpolate a profile defined by a lambda at the bottom of the domain.
 * It constantly (== Dirichlet) heats up the bottom of the domain and cools down towards the top. Also initializes
 * points that do not lie on the boundary.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp BC
 *
 * \section FullAppPlumeInCube-info Simulation info
 *
 * We add some code to print information of our setup like the number of DoFs and the coarse grid structure.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Info
 *
 * \section FullAppPlumeInCube-VTK VTK
 *
 * To visualize our results, we add the relevant function variables to our VTK output instance.
 * The output interval can directly be set in the constructor. Therefore it is not necessary
 * to count the number of time steps and check if VTK should be written. The respective member
 * is simply called in each time step and decides internally if output is written or not.
 * We also plot the initial conditions (time step == 0) and the coarse grid domain.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp VTK
 *
 * \section FullAppPlumeInCube-solvers Solvers
 *
 * To solve the Stokes equation we setup our geometric multigrid solver.
 * We need to define
 * - a coarse grid solver
 * - a smoother
 * - and the grid transfer operators.
 *
 * To easily setup a default case we can use the solver templates provied by HyTeG.
 * The components that are used are the following:
 *
 * For the coarse grid we employ a pressure preconditioned MinRes solver.
 *
 * To relax on the operator, we employ an inexact Uzawa iteration.
 * It requires to set an under-relaxation parameter for the pressure.
 * The optimal value can be calculated theoretically. However, we simply set it to
 * a fixed value wich appears to yield good convergence.
 *
 * The geometric multigrid solver has a parameter to specify an increased number of smoothing
 * steps on each coarser grid. From experiments, we find that this is required for the Uzawa-based
 * multigrid solver if a V-cycle is employed. We increase the number of pre- and post-smoothing
 * steps by 2 on each coarser level.
 *
 * The rest are standard multigrid components.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Solvers
 *
 * \section FullAppPlumeInCube-simulation Simulation loop
 *
 * Now we start the actual simulation. We define a short lambda to calculate the current residual
 * in the L2 norm.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Residual
 *
 * Then we start the main time stepping loop.
 *
 * First, we drive the Stokes flow by setting a force in z-direction that depends on the current temperature
 * field. Then we perform a few V-cycles to approximately solve the equation and calculate the average residual
 * reduction and some other infos.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Simulation Stokes
 *
 * With the result we advance the temperature transport by several time steps using the algebraic upwind operator.
 * We also write the VTK output.
 *
 * \snippet tutorials/06_FullAppPlumeInCube/Plume.cpp Simulation Advection
 *
 * This process is repeated until the simulation ends.
 *
 * Using the VTK data, the resulting plume rising up in the cube could be rendered like this:
 * \htmlonly
     <center>
     <table>
     <tr>
     <td align="center"><img src="plume.0000.png" width="80%"/></td>
     <td align="center"><img src="plume.0004.png" width="80%"/></td>
     </tr>
     <tr>
     <td align="center">initial state</td>
     <td align="center">after 400 time steps</td>
     </tr>
     <tr>
     <td align="center"><img src="plume.0007.png" width="80%"/></td>
     <td align="center"><img src="plume.0010.png" width="80%"/></td>
     </tr>
     <tr>
     <td align="center">after 700 time steps</td>
     <td align="center">after 1000 time steps</td>
     </tr>
     </table>
     </center>
     \endhtmlonly
 *
 *
 * \section FullAppPlumeInCube-fullApp Full Application
 * \include tutorials/06_FullAppPlumeInCube/Plume.cpp
 *
 */

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P1Transport.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   /// [Setup environment]

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   // Parameters

   auto config = std::make_shared< walberla::config::Config >();
   config->readParameterFile( "./PlumeParameters.prm" );

   const walberla::Config::BlockHandle mainConf = config->getBlock( "Parameters" );

   const uint_t minLevel   = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel   = mainConf.getParameter< uint_t >( "maxLevel" );
   const uint_t numVCycles = mainConf.getParameter< uint_t >( "numVCycles" );

   const real_t convectivity              = mainConf.getParameter< real_t >( "convectivity" );
   const real_t dt                        = mainConf.getParameter< real_t >( "dt" );
   const real_t plotDt                    = mainConf.getParameter< real_t >( "plotDt" );
   const real_t viscosity                 = mainConf.getParameter< real_t >( "viscosity" );
   const uint_t numStokesSteps            = mainConf.getParameter< uint_t >( "numStokesSteps" );
   const uint_t numTimeStepsPerStokesStep = mainConf.getParameter< uint_t >( "numTimeStepsPerStokesStep" );

   const bool vtkEnabled = mainConf.getParameter< bool >( "vtkEnabled" );

   /// [Setup environment]

   /// [Domain]

   MeshInfo              meshInfo = MeshInfo::meshCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   /// [Domain]

   /// [Discretization]

   P1P1StokesOperator       L( storage, minLevel, maxLevel );
   P1ConstantMassOperator M( storage, minLevel, maxLevel );
   P1Transport            transportOperator( storage, minLevel, maxLevel );

   P1StokesFunction< real_t > u( "u", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > f( "f", storage, minLevel, maxLevel );
   P1StokesFunction< real_t > r( "residual", storage, minLevel, maxLevel );
   P1Function< real_t >       temp( "temperature", storage, minLevel, maxLevel );

   /// [Discretization]

   /// [BC]

   std::function< real_t( const Point3D& ) > temperature = []( const Point3D& x ) {
      real_t temp_ = 1.0 - std::pow( x[2], 0.5 );
      return temp_ + 0.1 * ( 1.0 - x[2] ) * ( std::sin( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] ) );
   };

   temp.interpolate( temperature, maxLevel );

   /// [BC]

   /// [Info]

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo );

   const auto numDofsStokes      = numberOfGlobalDoFs< P1StokesFunctionTag >( *storage, maxLevel );
   const auto numDofsTemperature = numberOfGlobalDoFs< P1FunctionTag >( *storage, maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "Number of DoFs:" );
   WALBERLA_LOG_INFO_ON_ROOT( " - Stokes system (velocity + pressure): " << numDofsStokes );
   WALBERLA_LOG_INFO_ON_ROOT( " - Temperature:                         " << numDofsTemperature );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   /// [Info]

   /// [VTK]

   writeDomainPartitioningVTK( storage, "./vtk", "Plume_Domain" );

   VTKOutput vtkOutput( "./vtk", "StokesCubeTransport", storage, walberla::uint_c( std::ceil( plotDt / dt ) ) );

   vtkOutput.add( u );
   vtkOutput.add( temp );

   vtkOutput.write( maxLevel, 0 );

   /// [VTK]

   /// [Solvers]

   auto uzawaSolver = solvertemplates::stokesGMGUzawaSolver< P1P1StokesOperator >( storage, minLevel, maxLevel, 2, 2, 0.3 );

   /// [Solvers]

   // short lambda to calculate the residual of the discrete Stokes equation
   /// [Residual]
   auto calculateResidual = [&]() {
      L.apply( u, r, maxLevel, Inner );
      r.assign( { 1.0, -1.0 }, { f, r }, maxLevel, Inner );
      return sqrt( r.dotGlobal( r, maxLevel, Inner ) ) /
             real_c( numberOfGlobalDoFs< P1StokesFunctionTag >( *storage, maxLevel ) );
   };
   /// [Residual]

   real_t time                 = 0.0;
   real_t lastResidualL2       = 0;
   real_t avgResidualReduction = 0;

   // main simulation loop
   for ( uint_t stokesStep = 0; stokesStep < numStokesSteps; ++stokesStep )
   {
      /// [Simulation Stokes]
      // mass matrix applied to temperature and stored in right-hand side vector
      M.apply( temp, f.uvw()[2], maxLevel, All );
      f.uvw()[2].assign( { convectivity }, { f.uvw()[2] }, maxLevel, All );

      real_t currentResidualL2 = calculateResidual();

      // perform a number of V-cycles to approximately solve the Stokes system
      for ( uint_t i = 0; i < numVCycles; i++ )
      {
         uzawaSolver->solve( L, u, f, maxLevel );

         lastResidualL2    = currentResidualL2;
         currentResidualL2 = calculateResidual();
         avgResidualReduction += currentResidualL2 / lastResidualL2;
      }

      avgResidualReduction /= real_c( numVCycles );

      WALBERLA_LOG_INFO_ON_ROOT( "Simulated time: " << std::setw( 6 ) << std::setprecision( 2 ) << std::fixed << time
                                                    << "s (time step " << std::setw( 6 ) << stokesStep * numTimeStepsPerStokesStep
                                                    << ") | Stokes solver: final residual: " << std::scientific
                                                    << currentResidualL2 << " (avg. conv. rate: " << avgResidualReduction << ")" )
      /// [Simulation Stokes]

      /// [Simulation Advection]
      // advance the simulation a number of time steps with the current velocity solution
      for ( uint_t innerSteps = 0; innerSteps < numTimeStepsPerStokesStep; ++innerSteps )
      {
         time += dt;
         transportOperator.step( temp, u.uvw(), maxLevel, Inner, dt, viscosity );
         if ( vtkEnabled )
            vtkOutput.write( maxLevel, innerSteps + numTimeStepsPerStokesStep * stokesStep );
      }
      /// [Simulation Advection]
   }

   return EXIT_SUCCESS;
}
