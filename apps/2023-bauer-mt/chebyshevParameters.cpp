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

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/LaTeX/KeyValueStore.hpp"
#include "hyteg/dataexport/LaTeX/Table.hpp"
#include "hyteg/petsc/PETScManager.hpp"

#include "common.hpp"

using namespace hyteg;
using walberla::real_c;

void chebyshevParameters()
{
   const uint_t                  minLevel    = 2;
   const uint_t                  maxLevel    = 5;
   const std::array< real_t, 5 > lowerBounds = { 0.03, 0.05, 0.08, 0.12, 0.2 };
   const std::array< real_t, 5 > upperBounds = { 1.03, 1.05, 1.08, 1.12, 1.2 };

   const MeshInfo     cube = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
   const auto         zero = []( const Point3D& ) { return Point3D{ 0.0, 0.0, 0.0 }; };
   const n1e1::System system{
       cube,
       zero, // solution
       zero  // rhs
   };

   Params params{ "chebyshevParameters" };
   params.system                              = system;
   params.initialGuess                        = { []( const Point3D& ) {
      return Point3D{ real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                      real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                      real_c( walberla::math::realRandom( -1.0, 1.0 ) ) };
   } };
   params.maxLevel                            = maxLevel;
   params.computeAndStoreLocalElementMatrices = true;
   params.nMaxIterations                      = 12;

   latex::KeyValueStore store;
   params.store( store );
   WALBERLA_LOG_INFO_ON_ROOT( std::endl << store )
   store.writePgfKeys( "output", params.name );

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      params.maxLevel = level;
      latex::Table< 6 > table( { "lowerBound\\upperBound",
                                 walberla::format( "%f", upperBounds[0] ),
                                 walberla::format( "%f", upperBounds[1] ),
                                 walberla::format( "%f", upperBounds[2] ),
                                 walberla::format( "%f", upperBounds[3] ),
                                 walberla::format( "%f", upperBounds[4] ) } );

      for ( uint_t l = 0; l < lowerBounds.size(); ++l )
      {
         for ( uint_t u = 0; u < upperBounds.size(); ++u )
         {
            params.lowerBoundFactor = lowerBounds[l];
            params.upperBoundFactor = upperBounds[u];
            Results results         = solve( params );

            const real_t conv =
                std::pow( results.finalU2 / results.initU2, 1.0 / walberla::numeric_cast< real_t >( params.nMaxIterations ) );
            WALBERLA_LOG_INFO_ON_ROOT( "Level " << level << ": " << conv )

            table.addElement( l, 0, lowerBounds[l] );
            table.addElement( l, u + 1, conv );
         }
      }

      WALBERLA_LOG_INFO_ON_ROOT( std::endl << table )
      table.write( "output", walberla::format( "%s-level%i", params.name.c_str(), level ) );
   }
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   walberla::math::seedRandomGenerator( 0 );

   chebyshevParameters();
}
