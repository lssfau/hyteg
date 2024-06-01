/*
 * Copyright (c) 2024 Eugenio D'Ascoli.
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
#include "Convection.hpp"

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::WcTimer timeStepTimer;

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

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );
   terraneo::ConvectionSimulation      simulation( mainConf );
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
         return EXIT_SUCCESS;
      }
   }

   if ( simulation.getSimulationParams().simulationType == "CirculationModel" )
   {
      //if we reach this point in a circulation model, the max number of steps has been reached before the desired age
      WALBERLA_LOG_INFO_ON_ROOT( "Max timestep reached at age " << simulation.getSimulationParams().ageMa << " Ma." )

      WALBERLA_LOG_INFO_ON_ROOT( "Circulatin model ran from " << simulation.getSimulationParams().initialAge << " - "
                                                              << simulation.getSimulationParams().ageMa << " Ma." )
   }

   return EXIT_SUCCESS;
};