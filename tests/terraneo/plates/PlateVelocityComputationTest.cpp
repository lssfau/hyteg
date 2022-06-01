/*
 * Copyright (c) 2022 Berta Vilacis, Marcus Mohr.
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

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/extern/json.hpp"
#include "core/logging/Logging.h"

#include "hyteg/BuildInfo.hpp"
#include "hyteg/Git.hpp"

#include "terraneo/plates/PlateVelocityProvider.hpp"

using namespace hyteg;
using namespace terraneo;

typedef struct
{
   real_t age;
   real_t latitude;
   real_t longitude;
   uint_t plateID;
   vec3D  velocity;
} checkItem;

std::vector< checkItem > setupReferenceValues()
{
   std::vector< checkItem > refs;

   refs.push_back( {real_c( 1.0 ),
                    real_c( -31 ),
                    real_c( -56 ),
                    818,
                    {real_c( -8.269310829098373e-12 ), real_c( +2.604979937141433e-11 ), real_c( -1.383067664506167e-11 )}} );

   refs.push_back( {real_c( 1.0 ),
                    real_c( -35.0 ),
                    real_c( -56.0 ),
                    818,
                    {real_c( -6.110900513608632e-11 ), real_c( 1.782723509068086e-10 ), real_c( -1.027346527054028e-10 )}} );

   refs.push_back( {real_c( 23.0 ),
                    real_c( -90.0 ),
                    real_c( -25.0 ),
                    911,
                    {real_c( 2.123626519186899e-09 ), real_c( -5.857049309732643e-10 ), real_c( 1.256048277848174e-09 )}} );

   refs.push_back( {real_c( 34.0 ),
                    real_c( -30.0 ),
                    real_c( 75.0 ),
                    102,
                    {real_c( -5.79234735557856e-10 ), real_c( -1.21056872024404e-10 ), real_c( 1.18193342609513e-10 )}} );

   refs.push_back( {real_c( 50.0 ),
                    real_c( 10.0 ),
                    real_c( 20.0 ),
                    701,
                    {real_c( -4.363553948085126e-11 ), real_c( 3.345378170485724e-10 ), real_c( -4.153983837175462e-11 )}} );

   refs.push_back( {real_c( 82.0 ),
                    real_c( -150.0 ),
                    real_c( 0.0 ),
                    901,
                    {real_c( -8.292281411603950e-10 ), real_c( 1.436265271555702e-09 ), real_c( 2.479116222196206e-09 )}} );

   refs.push_back( {real_c( 83.0 ),
                    real_c( -60.0 ),
                    real_c( -25.0 ),
                    902,
                    {real_c( 1.765063140162903e-09 ), real_c( 3.200127619961332e-10 ), real_c( 1.298268167099986e-09 )}} );

   return refs;
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::string separator{"--------------------------------------------------"};
   WALBERLA_LOG_INFO_ON_ROOT( separator );
   printBuildInfo();
   printGitInfo();
   WALBERLA_LOG_INFO_ON_ROOT( separator );

   // initialise an oracle
   std::string                             dataDir{"../../data/terraneo/plates/"};
   std::string                             fnameTopologies      = dataDir + "topologies0-100Ma.geojson";
   std::string                             fnameReconstructions = dataDir + "Global_EarthByte_230-0Ma_GK07_AREPS.rot";
   terraneo::plates::PlateVelocityProvider oracle( fnameTopologies, fnameReconstructions );

   // set reference values for testing
   auto refVals = setupReferenceValues();

   // limits for deviation of computed values from reference
   real_t absTol = real_c( 1.0e-15 );
   real_t relTol = real_c( 3 * 1e-8 );

   // do the testing
   for ( uint_t idx = 0; idx < refVals.size(); ++idx )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " ======================================================\n"
                                 << "  Testing point with index " << idx << '\n'
                                 << " ----------------------------" );

      checkItem testCase{refVals[idx]};
      WALBERLA_LOG_INFO_ON_ROOT( " age = " << testCase.age );
      WALBERLA_LOG_INFO_ON_ROOT( " spherical coordinates = (" << testCase.latitude << ", " << testCase.longitude << " )" );
      std::vector< real_t > latLon{testCase.latitude, testCase.longitude};
      vec3D                 pointXYZ = terraneo::conversions::sph2cart( latLon );
      WALBERLA_LOG_INFO_ON_ROOT( " cartesian coordinates:\n" << pointXYZ );

      vec3D velXYZ = oracle.getPointVelocity( pointXYZ, testCase.age );
      WALBERLA_LOG_INFO_ON_ROOT( "" << velXYZ );

      // check component-wise differences
      for ( long k = 0; k < 3; ++k )
      {
         real_t absDiff = std::abs( testCase.velocity[k] - velXYZ[k] );
         real_t relDiff = absDiff / std::abs( testCase.velocity[k] );

         WALBERLA_LOG_INFO_ON_ROOT( "" << std::scientific << "\n computed value .... " << velXYZ[k] << "\n reference value ... "
                                       << testCase.velocity[k] << "\n absolute error .... " << absDiff
                                       << "\n relative .......... " << relDiff );

         WALBERLA_CHECK_LESS_EQUAL( absDiff, absTol );
         WALBERLA_CHECK_LESS_EQUAL( relDiff, relTol );
      }

      // orthogonality check
      real_t angleTol  = real_c( 1e-15 );
      real_t innerProd = velXYZ.dot( pointXYZ );
      real_t angle     = std::acos( innerProd / ( velXYZ.norm() * pointXYZ.norm() ) );
      angle *= real_c( 180 ) / terraneo::conversions::pi;
      real_t magDiff = std::abs( real_c( 90 ) - angle );

      WALBERLA_LOG_INFO_ON_ROOT( "" << std::scientific << "\n orthogonality check:"
                                 << "\n computed angle ...... " << angle
                                 << "\n error magnitude ..... " << magDiff );
      WALBERLA_CHECK_LESS_EQUAL( magDiff, angleTol );
   }

   return EXIT_SUCCESS;
}
