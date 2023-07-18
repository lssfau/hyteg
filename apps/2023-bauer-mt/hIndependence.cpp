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

#include "hyteg/dataexport/LaTeX/KeyValueStore.hpp"
#include "hyteg/dataexport/LaTeX/Table.hpp"
#include "hyteg/petsc/PETScManager.hpp"

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

   Params params{ "hIndependence" };
   params.maxLevel           = maxLevel;
   params.preSmoothSteps     = 3;
   params.postSmoothSteps    = 3;
   params.nMaxIterations     = 12;
   params.residual2Reduction = { 1e-6 };

   latex::KeyValueStore store;
   params.store( store );

   latex::Table< 13 > table( { "level",
                               "tet_poly_n_dofs",
                               "tet_poly_n_its",
                               "tet_poly_rate",
                               "tet_sine_n_dofs",
                               "tet_sine_n_its",
                               "tet_sine_rate",
                               "cube_poly_n_dofs",
                               "cube_poly_n_its",
                               "cube_poly_rate",
                               "cube_sine_n_dofs",
                               "cube_sine_n_its",
                               "cube_sine_rate" } );

   for ( uint_t i = 0; i < systems.size(); ++i )
   {
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         params.system   = systems[i];
         params.maxLevel = level;

         Results      results = solve( params );
         const real_t rate    = std::pow( results.finalResidual2 / results.initResidual2,
                                       1.0 / walberla::numeric_cast< real_t >( results.nIterations ) );
         WALBERLA_LOG_INFO_ON_ROOT( "Level " << level << std::endl
                                             << "  Number of DoFs            : " << results.numberOfGlobalDoFs << std::endl
                                             << "  Rate of residual reduction: " << rate << std::endl
                                             << "  Number of iterations      : " << results.nIterations )

         table.addElement( level - minLevel, 0, level );
         table.addElement( level - minLevel, 3 * i + 1, results.numberOfGlobalDoFs );
         table.addElement( level - minLevel, 3 * i + 2, results.nIterations );
         table.addElement( level - minLevel, 3 * i + 3, rate );
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
