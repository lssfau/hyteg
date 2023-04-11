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

#include "hyteg/petsc/PETScManager.hpp"

#include "KeyValueStore.hpp"
#include "Table.hpp"
#include "common.hpp"

using namespace hyteg;
using walberla::real_c;

void solverConvergenceCube()
{
   const uint_t                  minLevel = 3;
   const uint_t                  maxLevel = 6;
   const std::array< real_t, 3 > alphas   = { 0.01, 1.0, 100.0 };
   const std::array< real_t, 3 > betas    = { 0.01, 1.0, 100.0 };

   const MeshInfo     cube = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
   const auto         zero = []( const Point3D& ) { return Eigen::Vector3r{ 0.0, 0.0, 0.0 }; };
   const n1e1::System system{
       cube,
       zero, // solution
       zero  // rhs
   };

   Params params{ "solverConvergenceCube" };
   params.system                              = system;
   params.initialGuess                        = { []( const Point3D& ) {
      return Eigen::Vector3r{ real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                              real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                              real_c( walberla::math::realRandom( -1.0, 1.0 ) ) };
   } };
   params.maxLevel                            = maxLevel;
   params.computeAndStoreLocalElementMatrices = true;
   params.nMaxIterations                      = 12;

   KeyValueStore store;
   params.store( store );
   WALBERLA_LOG_INFO_ON_ROOT( std::endl << store )
   store.writePgfKeys( "output", params.name );

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      params.maxLevel = level;
      Table< 4 > table( { "alpha\\beta",
                          walberla::format( "%e", betas[0] ),
                          walberla::format( "%e", betas[1] ),
                          walberla::format( "%e", betas[2] ) } );

      for ( uint_t a = 0; a < alphas.size(); ++a )
      {
         for ( uint_t b = 0; b < betas.size(); ++b )
         {
            params.coefficients = { alphas[a], betas[b] };
            Results results     = solve( params );

            const real_t conv =
                std::pow( results.finalU2 / results.initU2, 1.0 / walberla::numeric_cast< real_t >( params.nMaxIterations ) );
            WALBERLA_LOG_INFO_ON_ROOT( "Level " << level << ": " << conv )

            table.addElement( a, 0, alphas[a] );
            table.addElement( a, b + 1, conv );
         }
      }

      WALBERLA_LOG_INFO_ON_ROOT( std::endl << table )
      table.write( "output", walberla::format( "%s-level%i", params.name.c_str(), level ) );
   }
}

void solverConvergenceTorus()
{
   const uint_t                  minLevel = 2;
   const uint_t                  maxLevel = 5;
   const real_t                  alpha    = 1.0;
   const std::array< real_t, 3 > betas    = { 0.01, 1.0, 100.0 };

   const MeshInfo     solidTorus = MeshInfo::meshTorus( 16, 8, 4.0, { 1.6 } );
   const auto         zero       = []( const Point3D& ) { return Eigen::Vector3r{ 0.0, 0.0, 0.0 }; };
   const n1e1::System system{
       solidTorus,
       zero, // solution
       zero  // rhs
   };

   Params params{ "solverConvergenceTorus" };
   params.coefficients                        = { alpha, 1.0 };
   params.system                              = system;
   params.initialGuess                        = { []( const Point3D& ) {
      return Eigen::Vector3r{ real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                              real_c( walberla::math::realRandom( -1.0, 1.0 ) ),
                              real_c( walberla::math::realRandom( -1.0, 1.0 ) ) };
   } };
   params.maxLevel                            = maxLevel;
   params.computeAndStoreLocalElementMatrices = true;
   params.nMaxIterations                      = 12;

   KeyValueStore store;
   params.store( store );
   WALBERLA_LOG_INFO_ON_ROOT( std::endl << store )
   store.writePgfKeys( "output", params.name );

   Table< 5 > table( { "level\\beta",
                       "n_dofs",
                       walberla::format( "%e", betas[0] ),
                       walberla::format( "%e", betas[1] ),
                       walberla::format( "%e", betas[2] ) } );

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      params.maxLevel = level;

      for ( uint_t b = 0; b < betas.size(); ++b )
      {
         params.coefficients = { alpha, betas[b] };
         Results results     = solve( params );

         const real_t conv =
             std::pow( results.finalU2 / results.initU2, 1.0 / walberla::numeric_cast< real_t >( params.nMaxIterations ) );
         WALBERLA_LOG_INFO_ON_ROOT( "Level " << level << ": " << conv )

         table.addElement( level - minLevel, 0, level );
         table.addElement( level - minLevel, 1, results.numberOfGlobalDoFs );
         table.addElement( level - minLevel, b + 2, conv );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( std::endl << table )
   table.write( "output", params.name );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   walberla::math::seedRandomGenerator( 0 );

   solverConvergenceCube();
   solverConvergenceTorus();
}
