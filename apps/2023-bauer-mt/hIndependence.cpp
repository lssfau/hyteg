/*
* Copyright (c) 2023 Daniel Bauer.
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

#include <vector>

#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "hyteg/petsc/PETScManager.hpp"

#include "KeyValueStore.hpp"
#include "Table.hpp"
#include "common.hpp"

using namespace hyteg;

void hIndependenceTest()
{
   const uint_t      minLevel = 3;
   const uint_t      maxLevel = 7;
   const std::vector systems  = { n1e1::System::polynomialOnSingleTet(),
                                  n1e1::System::sinusoidalOnSingleTet(),
                                  n1e1::System::polynomialOnCube(),
                                  n1e1::System::sinusoidalOnCube() };

   Params params{ "hIndependence",
                  // system
                  { 1.0, 1.0 }, // coefficients
                  systems[0],
                  {}, // initial guess = 0
                  // solver
                  0,        // min level
                  maxLevel, // max level
                  false,    // precompute element matrices
                  4,        // Chebyshev order
                  40,       // spectral estimation iterations
                  3,        // pre smooth
                  3,        // post smooth
                  10,       // max V-cycles
                  { 1e-6 }, // residual reduction
                  // output
                  false };

   KeyValueStore store;
   params.store( store );

   Table< 5 > table( { "level", "tet_poly", "tet_sine", "cube_poly", "cube_sine" } );

   for ( uint_t i = 0; i < systems.size(); ++i )
   {
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         params.system   = systems[i];
         params.maxLevel = level;

         Results results = solve( params );
         WALBERLA_LOG_INFO_ON_ROOT( "level " << level << ": " << results.nVCycles )

         table.addElement( level - minLevel, 0, level );
         table.addElement( level - minLevel, i + 1, results.nVCycles );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( std::endl << store )
   WALBERLA_LOG_INFO_ON_ROOT( std::endl << table )
   store.writePgfKeys( "output", params.name );
   table.write( "output", params.name );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   hIndependenceTest();
}
