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

#include <random>

#include "core/Environment.h"
#include "core/timing/TimingTree.h"

#include "forms/N1E1Form_mass_BInv_fp32.hpp"
#include "forms/N1E1Form_mass_BInv_fp64.hpp"
#include "forms/N1E1Form_mass_B_fp32.hpp"
#include "forms/N1E1Form_mass_B_fp64.hpp"
#include "forms/N1E1Form_mass_LU_fp32.hpp"
#include "forms/N1E1Form_mass_LU_fp64.hpp"

using namespace hyteg;

/// returns a random point on the unit sphere, uniformally wrt the angle
PointND< float, 3 > rndOnUnitSphere()
{
   static std::mt19937                      generator;
   static std::normal_distribution< float > normalDist{ 0.0f };

   const float x{ normalDist( generator ) };
   const float y{ normalDist( generator ) };
   const float z{ normalDist( generator ) };

   const float rInv = 1.0f / std::sqrt( x * x + y * y + z * z );

   return { { x * rInv, y * rInv, z * rInv } };
}

void testUnitSphere()
{
   const n1e1::N1E1Form_mass_BInv_fp32 formBInv32;
   const n1e1::N1E1Form_mass_BInv_fp64 formBInv64;
   const n1e1::N1E1Form_mass_B_fp32    formB32;
   const n1e1::N1E1Form_mass_B_fp64    formB64;
   const n1e1::N1E1Form_mass_LU_fp32   formLU32;
   const n1e1::N1E1Form_mass_LU_fp64   formLU64;

   Matrix< float, 6, 6 >  matBInv32;
   Matrix< double, 6, 6 > matBInv64;
   Matrix< float, 6, 6 >  matB32;
   Matrix< double, 6, 6 > matB64;
   Matrix< float, 6, 6 >  matLU32;
   Matrix< double, 6, 6 > matLU64;

   double frobeniusBInv{ 0.0 };
   double frobeniusB{ 0.0 };
   double frobeniusLU{ 0.0 };

   double maxBInv{ 0.0 };
   double maxB{ 0.0 };
   double maxLU{ 0.0 };

   const int n = 1000;
   for ( int cell = 0; cell < n; ++cell )
   {
      const std::array< PointND< float, 3 >, 4 > coords32{
          PointND< float, 3 >{ { 0, 0, 0 } }, rndOnUnitSphere(), rndOnUnitSphere(), rndOnUnitSphere() };
      std::array< PointND< double, 3 >, 4 > coords64;
      for ( uint_t i = 0; i < coords32.size(); ++i )
      {
         coords64[i].vector_ = coords32[i].vector_.cast< double >();
      }

      const std::array< walberla::int16_t, 6 > edgeDirections{ 1, 1, 1, 1, 1, 1 };

      formBInv32.integrateAll( coords32, edgeDirections, matBInv32 );
      formBInv64.integrateAll( coords64, edgeDirections, matBInv64 );
      formB32.integrateAll( coords32, edgeDirections, matB32 );
      formB64.integrateAll( coords64, edgeDirections, matB64 );
      formLU32.integrateAll( coords32, edgeDirections, matLU32 );
      formLU64.integrateAll( coords64, edgeDirections, matLU64 );

      frobeniusBInv += ( matBInv64 - matBInv32.cast< double >() ).norm();
      frobeniusB += ( matB64 - matB32.cast< double >() ).norm();
      frobeniusLU += ( matLU64 - matLU32.cast< double >() ).norm();
      maxBInv = std::max( maxBInv, ( matBInv64 - matBInv32.cast< double >() ).norm() );
      maxB    = std::max( maxB, ( matB64 - matB32.cast< double >() ).norm() );
      maxLU   = std::max( maxLU, ( matLU64 - matLU32.cast< double >() ).norm() );
   }

   frobeniusBInv /= n;
   frobeniusB /= n;
   frobeniusLU /= n;

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( frobeniusBInv )
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( frobeniusB )
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( frobeniusLU )

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( maxBInv )
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( maxB )
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT( maxLU )
}

void testTime()
{
   const n1e1::N1E1Form_mass_BInv_fp64 formBInv64;
   const n1e1::N1E1Form_mass_B_fp64    formB64;
   const n1e1::N1E1Form_mass_LU_fp64   formLU64;

   Matrix< double, 6, 6 > matBInv64;
   Matrix< double, 6, 6 > matB64;
   Matrix< double, 6, 6 > matLU64;

   Eigen::Matrix< double, 3, 4 > maxCoords;

   const std::array< PointND< float, 3 >, 4 > coords32{
       PointND< float, 3 >{ { 0, 0, 0 } }, rndOnUnitSphere(), rndOnUnitSphere(), rndOnUnitSphere() };
   std::array< PointND< double, 3 >, 4 > coords64;
   for ( uint_t i = 0; i < coords32.size(); ++i )
   {
      coords64[i].vector_ = coords32[i].vector_.cast< double >();
   }

   const std::array< walberla::int16_t, 6 > edgeDirections{ 1, 1, 1, 1, 1, 1 };
   const int                                n = 1e6;
   walberla::WcTimingTree                   tt;

   tt.start( "BInv64" );
   for ( int cell = 0; cell < n; ++cell )
   {
      formBInv64.integrateAll( coords64, edgeDirections, matBInv64 );
   }
   tt.stop( "BInv64" );

   tt.start( "B64" );
   for ( int cell = 0; cell < n; ++cell )
   {
      formB64.integrateAll( coords64, edgeDirections, matB64 );
   }
   tt.stop( "B64" );

   tt.start( "LU64" );
   for ( int cell = 0; cell < n; ++cell )
   {
      formLU64.integrateAll( coords64, edgeDirections, matLU64 );
   }
   tt.stop( "LU64" );

   WALBERLA_LOG_INFO_ON_ROOT( tt )
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );

   testUnitSphere();
   testTime();
}
