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

#include "terraneo/dataimport/ParameterIO.hpp"

/**
 * @file TerraNeoParameterTest.cpp
 * @brief Test program for data input/output modularization.
 *
 * The TerraNeoParameterTest program is used to test the functionality of the TerraNeo data input/output module.
 * It reads a parameter file, parses the configuration, and prints the configuration settings.
 *
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line argument strings.
 * @return 0 on successful execution.
 */

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "../../data/param/TerraNeoParameterTest.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );
   // If a viscosity file is given, parse config will implicitly call readDataFile and
   // will store the viscosity in the corresponding data structure
   terraneo::parseConfig( mainConf );

   bool loggingToFile = true;

   if ( loggingToFile )
   {
      terraneo::printConfig();
   }

   if ( terraneo::simulationParam.haveViscosityProfile )
   {
      // Check entries in radius column of profile file (non-dimensional values)

      WALBERLA_CHECK_FLOAT_EQUAL( terraneo::physicalParam.viscosityProfile[0][0], terraneo::domainParam.rMax );
      WALBERLA_CHECK_FLOAT_EQUAL( terraneo::physicalParam.viscosityProfile.back()[0], terraneo::domainParam.rMin );

      // check if second data column is non-zero
      for ( size_t i = 0; i < terraneo::physicalParam.viscosityProfile.size(); i++ )
      {
         WALBERLA_CHECK_GREATER( terraneo::physicalParam.viscosityProfile[i][1], terraneo::real_c( 0 ) );
      }
   }
   return 0;
}