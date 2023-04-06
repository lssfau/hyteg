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

#include <cmath>

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/Environment.h"

#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/mesh/MeshInfo.hpp"

#include "KeyValueStore.hpp"
#include "Table.hpp"
#include "common.hpp"

using namespace hyteg;
using walberla::real_c;

void smoothingFactor( const uint_t n1e1SmoothSteps,
                      const uint_t p1SmoothSteps,
                      const real_t errorReductionFactor,
                      const uint_t nMaxIts )
{
   using walberla::math::pi;

   const uint_t level = 4;
   const uint_t minK  = 1;
   const uint_t maxK  = 16;

   std::string smootherDesc = walberla::format( "%i_%i", n1e1SmoothSteps, p1SmoothSteps );

   const MeshInfo     solidTorus = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );
   const auto         zero       = []( const Point3D& ) { return Eigen::Vector3r{ 0.0, 0.0, 0.0 }; };
   const n1e1::System system{
       solidTorus,
       zero, // solution
       zero  // rhs
   };

   Params params{ "smoothingFactor/" + smootherDesc };
   params.system                              = system;
   params.minLevel                            = level;
   params.maxLevel                            = level;
   params.computeAndStoreLocalElementMatrices = true;
   params.n1e1SmoothSteps                     = n1e1SmoothSteps;
   params.p1SmoothSteps                       = p1SmoothSteps;
   params.nMaxIterations                      = nMaxIts;
   params.u2Reduction                         = { 1.0 / errorReductionFactor };

   KeyValueStore store;
   params.store( store );

   Table< 4 > table( { "k", smootherDesc + "_fourier", smootherDesc + "_curl-free", smootherDesc + "_div-free" } );

   auto fourierMode = []( const int k, const Point3D& p ) {
      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];

      const real_t val = std::sin( real_c( k ) * x * pi ) * std::sin( real_c( k ) * y * pi ) * std::sin( real_c( k ) * z * pi );
      return Eigen::Vector3r{ val, val, val };
   };

   auto curlFreeMode = []( const int k, const Point3D& p ) {
      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];

      return Eigen::Vector3r{
          std::cos( real_c( k ) * x * pi ) * std::sin( real_c( k ) * y * pi ) * std::sin( real_c( k ) * z * pi ),
          std::sin( real_c( k ) * x * pi ) * std::cos( real_c( k ) * y * pi ) * std::sin( real_c( k ) * z * pi ),
          std::sin( real_c( k ) * x * pi ) * std::sin( real_c( k ) * y * pi ) * std::cos( real_c( k ) * z * pi ) };
   };

   auto divFreeMode = []( const int k, const Point3D& p ) {
      using std::sin;
      using std::cos;
      const real_t kp = real_c( k ) * pi;

      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];

      return Eigen::Vector3r{ sin( x * kp ) * cos( y * kp ) * sin( z * kp ) - sin( x * kp ) * sin( y * kp ) * cos( z * kp ),
                              sin( x * kp ) * sin( y * kp ) * cos( z * kp ) - cos( x * kp ) * sin( y * kp ) * sin( z * kp ),
                              cos( x * kp ) * sin( y * kp ) * sin( z * kp ) - sin( x * kp ) * cos( y * kp ) * sin( z * kp ) };
   };

   WALBERLA_LOG_INFO_ON_ROOT( "╭──────────────────────────────────╮" )
   WALBERLA_LOG_INFO_ON_ROOT( "│ Smoothing steps in 𝒩 (curl)^⊥: " << n1e1SmoothSteps << " │" );
   WALBERLA_LOG_INFO_ON_ROOT( "│ Smoothing steps in 𝒩 (curl)  : " << p1SmoothSteps << " │" );
   WALBERLA_LOG_INFO_ON_ROOT( "╰──────────────────────────────────╯" )

   WALBERLA_LOG_INFO_ON_ROOT( "Fourier mode" )
   for ( uint_t k = minK; k < maxK; ++k )
   {
      params.initialGuess = std::bind( fourierMode, k, std::placeholders::_1 );
      Results results     = smooth( params );

      WALBERLA_LOG_INFO_ON_ROOT( "  k = " << std::setw( 2 ) << k << ": " << std::setw( 3 ) << results.nIterations );

      table.addElement( k - minK, 0, k );
      table.addElement( k - minK, 1, results.nIterations );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Curl-free mode" )
   for ( uint_t k = minK; k < maxK; ++k )
   {
      params.initialGuess = std::bind( curlFreeMode, k, std::placeholders::_1 );
      Results results     = smooth( params );

      WALBERLA_LOG_INFO_ON_ROOT( "  k = " << std::setw( 2 ) << k << ": " << std::setw( 3 ) << results.nIterations );

      table.addElement( k - minK, 0, k );
      table.addElement( k - minK, 2, results.nIterations );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Div-free mode" )
   for ( uint_t k = minK; k < maxK; ++k )
   {
      params.initialGuess = std::bind( divFreeMode, k, std::placeholders::_1 );
      Results results     = smooth( params );

      WALBERLA_LOG_INFO_ON_ROOT( "  k = " << std::setw( 2 ) << k << ": " << std::setw( 3 ) << results.nIterations );

      table.addElement( k - minK, 0, k );
      table.addElement( k - minK, 3, results.nIterations );
   }

   WALBERLA_LOG_INFO_ON_ROOT( std::endl << store )
   WALBERLA_LOG_INFO_ON_ROOT( std::endl << table )

   store.writePgfKeys( "output", "smoothingFactor-" + smootherDesc );
   table.write( "output", "smoothingFactor-" + smootherDesc );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   smoothingFactor( 1, 0, 5.0, 100 );
   smoothingFactor( 1, 1, 100.0, 100 );
}
