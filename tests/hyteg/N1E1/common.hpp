/*
* Copyright (c) 2022-2023 Daniel Bauer.
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

#pragma once

// This file contains some common functions and definitions used in the N1E1 tests

#include "core/DataTypes.h"
#include "core/math/Constants.h"

#include "hyteg/geometry/TokamakMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {
namespace n1e1 {

using VectorField = std::function< Point3D( const Point3D& ) >;

/// Some example problems (domain, analytical solution, right hand side) for
///   curl curl u + u = f  in Ω
///             u × n = 0  on ∂Ω
class System
{
 public:
   MeshInfo                                        domain_;
   VectorField                                     analyticalSol_;
   VectorField                                     rhs_;
   std::function< void( SetupPrimitiveStorage& ) > setMap_;

   System( MeshInfo domain, VectorField analyticalSol, VectorField rhs )
   : domain_( domain )
   , analyticalSol_( analyticalSol )
   , rhs_( rhs )
   , setMap_( []( SetupPrimitiveStorage& ) {} )
   {}

   System( MeshInfo domain, VectorField analyticalSol, VectorField rhs, std::function< void( SetupPrimitiveStorage& ) > setMap )
   : domain_( domain )
   , analyticalSol_( analyticalSol )
   , rhs_( rhs )
   , setMap_( setMap )
   {}

   static System polynomialOnSingleTet()
   {
      return System{
          MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/tet_1el.msh" ) ),

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ y * z * ( x + y + z - 1 ), x * z * ( x + y + z - 1 ), x * y * ( x + y + z - 1 ) };
          },

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{
                 y * z * ( x + y + z - 1 ) - y - z, x * z * ( x + y + z - 1 ) - x - z, x * y * ( x + y + z - 1 ) - x - y };
          },
      };
   }

   static System sinusoidalOnSingleTet()
   {
      using std::sin;
      using walberla::real_c;

      return System{
          MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/tet_1el_variant.msh" ) ),

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ sin( y ) * sin( z ) * sin( x + y + z - 1 ),
                             sin( x ) * sin( z ) * sin( x + y + z - 1 ),
                             sin( x ) * sin( y ) * sin( x + y + z - 1 ) };
          },

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ real_c( -sin( x ) * sin( y ) * sin( x + y + z - 1 ) - sin( x ) * sin( z ) * sin( x + y + z - 1 ) +
                                     5 * sin( y ) * sin( z ) * sin( x + y + z - 1 ) + sin( y ) * cos( x ) * cos( x + y + z - 1 ) -
                                     2 * sin( y ) * cos( z ) * cos( x + y + z - 1 ) + sin( z ) * cos( x ) * cos( x + y + z - 1 ) -
                                     2 * sin( z ) * cos( y ) * cos( x + y + z - 1 ) ),
                             real_c( -sin( x ) * sin( y ) * sin( x + y + z - 1 ) +
                                     5 * sin( x ) * sin( z ) * sin( x + y + z - 1 ) + sin( x ) * cos( y ) * cos( x + y + z - 1 ) -
                                     2 * sin( x ) * cos( z ) * cos( x + y + z - 1 ) - sin( y ) * sin( z ) * sin( x + y + z - 1 ) -
                                     2 * sin( z ) * cos( x ) * cos( x + y + z - 1 ) +
                                     sin( z ) * cos( y ) * cos( x + y + z - 1 ) ),
                             real_c( 5 * sin( x ) * sin( y ) * sin( x + y + z - 1 ) - sin( x ) * sin( z ) * sin( x + y + z - 1 ) -
                                     2 * sin( x ) * cos( y ) * cos( x + y + z - 1 ) + sin( x ) * cos( z ) * cos( x + y + z - 1 ) -
                                     sin( y ) * sin( z ) * sin( x + y + z - 1 ) - 2 * sin( y ) * cos( x ) * cos( x + y + z - 1 ) +
                                     sin( y ) * cos( z ) * cos( x + y + z - 1 ) ) };
          },
      };
   }

   static System polynomialOnPyramid()
   {
      using walberla::real_c;

      return System{
          MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/pyramid_2el_variant.msh" ) ),

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ z * ( z - 2 * x ) * ( z + 2 * x - 2 ) * ( z - 2 * y ) * ( z + 2 * y - 2 ),
                             z * ( z - 2 * x ) * ( z + 2 * x - 2 ) * ( z - 2 * y ) * ( z + 2 * y - 2 ),
                             1 * ( z - 2 * x ) * ( z + 2 * x - 2 ) * ( z - 2 * y ) * ( z + 2 * y - 2 ) };
          },

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{
                 real_c( 16 * x * x * y * y * z - 16 * x * x * y * z - 4 * x * x * std::pow( z, 3 ) + 8 * x * x * z * z -
                         8 * x * x * z - 16 * x * x - 16 * x * y * y * z + 80 * x * y * z + 4 * x * std::pow( z, 3 ) -
                         8 * x * z * z - 40 * x * z + 32 * x - 4 * y * y * std::pow( z, 3 ) + 8 * y * y * z * z + 24 * y * y * z -
                         16 * y * y + 4 * y * std::pow( z, 3 ) - 8 * y * z * z - 56 * y * z + 16 * y + std::pow( z, 5 ) -
                         4 * std::pow( z, 4 ) - 8 * std::pow( z, 3 ) + 32 * z * z - 8 ),
                 real_c( 16 * x * x * y * y * z - 16 * x * x * y * z - 4 * x * x * std::pow( z, 3 ) + 8 * x * x * z * z +
                         24 * x * x * z - 16 * x * x - 16 * x * y * y * z + 80 * x * y * z + 4 * x * std::pow( z, 3 ) -
                         8 * x * z * z - 56 * x * z + 16 * x - 4 * y * y * std::pow( z, 3 ) + 8 * y * y * z * z - 8 * y * y * z -
                         16 * y * y + 4 * y * std::pow( z, 3 ) - 8 * y * z * z - 40 * y * z + 32 * y + std::pow( z, 5 ) -
                         4 * std::pow( z, 4 ) - 8 * std::pow( z, 3 ) + 32 * z * z - 8 ),
                 real_c( 16 * x * x * y * y + 16 * x * x * y - 4 * x * x * z * z + 8 * x * x * z - 48 * x * x + 16 * x * y * y -
                         48 * x * y - 20 * x * z * z + 24 * x * z + 48 * x - 4 * y * y * z * z + 8 * y * y * z - 48 * y * y -
                         20 * y * z * z + 24 * y * z + 48 * y + std::pow( z, 4 ) - 4 * std::pow( z, 3 ) + 44 * z * z - 64 * z ) };
          },
      };
   }

   static System polynomialOnCube()
   {
      return System{
          MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 ),

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ y * ( 1 - y ) * z * ( 1 - z ), x * ( 1 - x ) * z * ( 1 - z ), x * ( 1 - x ) * y * ( 1 - y ) };
          },

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ 2 * ( y * ( 1 - y ) + z * ( 1 - z ) ) + y * ( 1 - y ) * z * ( 1 - z ),
                             2 * ( x * ( 1 - x ) + z * ( 1 - z ) ) + x * ( 1 - x ) * z * ( 1 - z ),
                             2 * ( x * ( 1 - x ) + y * ( 1 - y ) ) + x * ( 1 - x ) * y * ( 1 - y ) };
          },
      };
   }

   static System sinusoidalOnCube()
   {
      using std::sin;
      using walberla::math::pi;

      return System{
          MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), 1, 1, 1 ),

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ sin( 2 * pi * y ) * sin( 2 * pi * z ),
                             sin( 2 * pi * x ) * sin( 2 * pi * z ),
                             sin( 2 * pi * x ) * sin( 2 * pi * y ) };
          },

          []( const Point3D& p ) {
             const real_t x = p[0];
             const real_t y = p[1];
             const real_t z = p[2];
             return Point3D{ sin( 2 * pi * y ) * sin( 2 * pi * z ) + 8 * pi * pi * sin( 2 * pi * y ) * sin( 2 * pi * z ),
                             sin( 2 * pi * x ) * sin( 2 * pi * z ) + 8 * pi * pi * sin( 2 * pi * x ) * sin( 2 * pi * z ),
                             sin( 2 * pi * x ) * sin( 2 * pi * y ) + 8 * pi * pi * sin( 2 * pi * x ) * sin( 2 * pi * y ) };
          },
      };
   }

   static System onToroidalSlice()
   {
      const uint_t                toroidalResolution         = 34;
      const uint_t                poloidalResolution         = 6;
      const real_t                radiusOriginToCenterOfTube = 2;
      const std::vector< real_t > tubeLayerRadii             = { 0.4 };
      const real_t                torodialStartAngle         = 0.0;
      const real_t                polodialStartAngle         = 0.0;
      const uint_t                numToroidalSlices          = 1;
      const real_t                delta                      = 0;
      const real_t                r1                         = tubeLayerRadii.back();
      const real_t                r2                         = tubeLayerRadii.back();

      const real_t R = radiusOriginToCenterOfTube;
      const real_t r = tubeLayerRadii.back();

      return System{
          MeshInfo::meshTorus( toroidalResolution,
                               poloidalResolution,
                               radiusOriginToCenterOfTube,
                               tubeLayerRadii,
                               torodialStartAngle,
                               polodialStartAngle,
                               numToroidalSlices ),

          [R, r]( const Point3D& xVec ) {
             const real_t x    = xVec[0];
             const real_t y    = xVec[1];
             const real_t z    = xVec[2];
             const real_t tmp0 = std::sqrt( std::pow( x, 2 ) + std::pow( y, 2 ) );
             const real_t tmp1 = ( -std::pow( r, 2 ) + std::pow( z, 2 ) + std::pow( -R + tmp0, 2 ) ) / tmp0;
             const real_t u0   = tmp1 * y;
             const real_t u1   = -tmp1 * x;
             const real_t u2   = 0;
             return Point3D{ u0, u1, u2 };
          },

          [R, r]( const Point3D& xVec ) {
             const real_t x     = xVec[0];
             const real_t y     = xVec[1];
             const real_t z     = xVec[2];
             const real_t tmp0  = std::pow( x, 2 );
             const real_t tmp1  = std::pow( y, 2 );
             const real_t tmp2  = tmp0 + tmp1;
             const real_t tmp3  = std::sqrt( tmp2 );
             const real_t tmp4  = 1.0 / tmp3;
             const real_t tmp5  = std::pow( y, 3 );
             const real_t tmp6  = std::pow( tmp2, -3.0 / 2.0 );
             const real_t tmp7  = 2 * tmp6;
             const real_t tmp8  = tmp0 * y;
             const real_t tmp9  = -R + tmp3;
             const real_t tmp10 = 8 / tmp2;
             const real_t tmp11 = std::pow( tmp2, -2 );
             const real_t tmp12 = -std::pow( r, 2 ) + std::pow( tmp9, 2 ) + std::pow( z, 2 );
             const real_t tmp13 = 3 * tmp12 / std::pow( tmp2, 5.0 / 2.0 );
             const real_t tmp14 = tmp4 * x;
             const real_t tmp15 = std::pow( x, 3 );
             const real_t tmp16 = tmp1 * x;
             const real_t tmp17 = 6 * tmp11 * tmp9;
             const real_t u0    = 6 * tmp0 * tmp11 * tmp9 * y - tmp10 * tmp9 * y + 6 * tmp11 * tmp5 * tmp9 + tmp12 * tmp4 * y +
                               4 * tmp12 * tmp6 * y - tmp13 * tmp5 - tmp13 * tmp8 - 2 * tmp4 * y - tmp5 * tmp7 - tmp7 * tmp8;
             const real_t u1 = tmp10 * tmp9 * x - tmp12 * tmp14 - 4 * tmp12 * tmp6 * x + tmp13 * tmp15 + tmp13 * tmp16 +
                               2 * tmp14 - tmp15 * tmp17 + tmp15 * tmp7 - tmp16 * tmp17 + tmp16 * tmp7;
             const real_t u2 = 0;
             return Point3D{ u0, u1, u2 };
          },

          [=]( SetupPrimitiveStorage& setupStorage ) {
             TokamakMap::setMap( setupStorage,
                                 toroidalResolution,
                                 poloidalResolution,
                                 radiusOriginToCenterOfTube,
                                 tubeLayerRadii,
                                 torodialStartAngle,
                                 polodialStartAngle,
                                 delta,
                                 r1,
                                 r2 );
          },
      };
   }
};

/// Test that the L2 convergence factor is ≤ 1/4 (+1/40).
void L2ConvergenceTest( const uint_t                                                                             minLevel,
                        const uint_t                                                                             maxLevel,
                        const System                                                                             system,
                        std::function< real_t( const uint_t level, const System& system, const bool writeVtk ) > test,
                        const bool writeVTK = false );

} // namespace n1e1
} // namespace hyteg
