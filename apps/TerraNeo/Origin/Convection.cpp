/*
 * Copyright (c) 2024 Eugenio D'Ascoli, Andreas Burkhart,
 * Nils Kohl, Hamish Brown, Ponsuganth Ilangovan, Marcus Mohr.
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

#include "hyteg/Git.hpp"

#include "Simulation.hpp"

template < typename ConvectionSimulation_T >
void startSimulation( ConvectionSimulation_T simulation )
{
   walberla::WcTimer timeStepTimer;
   hyteg::buildinfo::printGitInfo();

   simulation.init();

   // Starting simulation

   for ( uint_t i = 0; i < simulation.getSimulationParams().maxNumTimesteps; ++i )
   {
      timeStepTimer.start();
      simulation.step();
      timeStepTimer.end();
      real_t timeStepTotal = timeStepTimer.last();
      WALBERLA_LOG_INFO_ON_ROOT( "Total time for step: " << timeStepTotal )

      if ( simulation.getSimulationParams().simulationType == "CirculationModel" &&
           ( simulation.getSimulationParams().ageMa <= simulation.getSimulationParams().finalAge ) )
      {
         //triggered following the early return from convSim.step(), due to current timestep (ageMa) passing the user-defined final time
         WALBERLA_LOG_INFO_ON_ROOT( "Circulation model finished at " << simulation.getSimulationParams().agePrev << " Ma." )
         return;
      }
      // Check if the simulation has reached the desired age for a convection model and finish the simulation
      if ( simulation.getSimulationParams().simulationType == "ConvectionModel" &&
           ( simulation.getSimulationParams().modelRunTimeMa >= simulation.getSimulationParams().finalAge ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Convection model finished at " << simulation.getSimulationParams().modelRunTimeMa
                                                                    << " Ma. After: " << simulation.getSimulationParams().timeStep
                                                                    << " timesteps." )
         return;
      }
   }

   if ( simulation.getSimulationParams().simulationType == "CirculationModel" )
   {
      //if we reach this point in a circulation model, the max number of steps has been reached before the desired age
      WALBERLA_LOG_INFO_ON_ROOT( "Max timestep reached at age " << simulation.getSimulationParams().ageMa << " Ma." )

      WALBERLA_LOG_INFO_ON_ROOT( "Circulation model ran from " << simulation.getSimulationParams().initialAge << " - "
                                                               << simulation.getSimulationParams().ageMa << " Ma." )
   }
   if ( simulation.getSimulationParams().timingAnalysis )
   {
      simulation.outputTimingTree();
   }
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager manager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./parameters.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle feSpacesConf = cfg->getBlock( "FunctionSpaces" );
   const walberla::Config::BlockHandle mainConf     = cfg->getBlock( "Parameters" );

   std::string temperatureSpace = feSpacesConf.getParameter< std::string >( "TemperatureSpace" );
   std::string viscositySpace   = feSpacesConf.getParameter< std::string >( "ViscositySpace" );

   if ( temperatureSpace == "P2" && viscositySpace == "P0" )
   {
      using TemperatureFunction_T = P2Function< real_t >;
      using ViscosityFunction_T   = P0Function< real_t >;

      terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > simulation( mainConf );
      startSimulation< terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > >( simulation );
   }
   else if ( temperatureSpace == "P2" && viscositySpace == "P1" )
   {
      using TemperatureFunction_T = P2Function< real_t >;
      using ViscosityFunction_T   = P1Function< real_t >;

      terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > simulation( mainConf );
      startSimulation< terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > >( simulation );
   }
   else if ( temperatureSpace == "P1" && viscositySpace == "P0" )
   {
      using TemperatureFunction_T = P1Function< real_t >;
      using ViscosityFunction_T   = P0Function< real_t >;

      terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > simulation( mainConf );
      startSimulation< terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > >( simulation );
   }
   else if ( temperatureSpace == "P1" && viscositySpace == "P1" )
   {
      using TemperatureFunction_T = P1Function< real_t >;
      using ViscosityFunction_T   = P1Function< real_t >;

      terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > simulation( mainConf );
      startSimulation< terraneo::ConvectionSimulation< TemperatureFunction_T, ViscosityFunction_T > >( simulation );
   }
   else
   {
      WALBERLA_ABORT( "This combination of function spaces is not available" );
   }

   return EXIT_SUCCESS;
}
