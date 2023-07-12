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

#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/KeyValueStore.hpp"
#include "hyteg/dataexport/Table.hpp"
#include "hyteg/petsc/PETScManager.hpp"

#include "common.hpp"

using namespace hyteg;
using walberla::real_t;

void L2ConvergenceTest()
{
   const uint_t      minLevel = 4;
   const uint_t      maxLevel = 7;
   const std::vector systems  = { n1e1::System::polynomialOnSingleTet(),
                                  n1e1::System::sinusoidalOnSingleTet(),
                                  n1e1::System::polynomialOnCube(),
                                  n1e1::System::sinusoidalOnCube() };

   Params params{ "L2Convergence" };
   params.maxLevel           = maxLevel;
   params.preSmoothSteps     = 3;
   params.postSmoothSteps    = 3;
   params.nMaxIterations     = 40;
   params.residual2Reduction = { 1e-11 };

   KeyValueStore store;
   params.store( store );

   Table< 5 > l2Error( { "level", "tet_poly", "tet_sine", "cube_poly", "cube_sine" } );
   Table< 5 > convergenceFactor( { "level", "tet_poly", "tet_sine", "cube_poly", "cube_sine" } );

   for ( uint_t i = 0; i < systems.size(); ++i )
   {
      params.system   = systems[i];
      params.maxLevel = minLevel;
      Results results = solve( params );
      real_t  err     = results.finalErrL2;

      WALBERLA_LOG_INFO_ON_ROOT( "L2-error on level " << minLevel << ": " << std::scientific << err );

      l2Error.addElement( 0, 0, minLevel );
      l2Error.addElement( 0, i + 1, err );

      for ( uint_t level = minLevel + 1; level <= maxLevel; level++ )
      {
         params.maxLevel           = level;
         Results      resultsFiner = solve( params );
         const real_t errFiner     = resultsFiner.finalErrL2;
         const real_t computedRate = errFiner / err;

         WALBERLA_LOG_INFO_ON_ROOT( "L2-error on level " << level << ": " << std::scientific << errFiner );
         WALBERLA_LOG_INFO_ON_ROOT( "Convergence rate level " << level << "/" << level - 1 << ": " << computedRate );

         l2Error.addElement( level - minLevel, 0, level );
         l2Error.addElement( level - minLevel, i + 1, errFiner );

         convergenceFactor.addElement( level - minLevel - 1, 0, walberla::format( "%i/%i", level, level - 1 ) );
         convergenceFactor.addElement( level - minLevel - 1, i + 1, computedRate );

         err = errFiner;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( std::endl << store )
   WALBERLA_LOG_INFO_ON_ROOT( std::endl << l2Error )
   WALBERLA_LOG_INFO_ON_ROOT( std::endl << convergenceFactor )
   store.writePgfKeys( "output", params.name );
   l2Error.write( "output", params.name + "-L2Error" );
   convergenceFactor.write( "output", params.name + "-convergenceFactor" );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   L2ConvergenceTest();
}
