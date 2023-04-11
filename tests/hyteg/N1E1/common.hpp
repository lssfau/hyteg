/*
* Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/mesh/MeshInfo.hpp"

namespace hyteg {
namespace n1e1 {

using VectorField = std::function< Point3D( const Point3D& ) >;

/// Some example problems (domain, analytical solution, right hand side) for
///   curl curl u + u = f  in Ω
///             u × n = 0  on ∂Ω
class System
{
 public:
   MeshInfo    domain_;
   VectorField analyticalSol_;
   VectorField rhs_;

   System( MeshInfo domain, VectorField analyticalSol, VectorField rhs )
   : domain_( domain )
   , analyticalSol_( analyticalSol )
   , rhs_( rhs )
   {}

   static System polynomialOnSingleTet()
   {
      return System{
          MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ),

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
          MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ),

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
          MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),

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
};

/// Test that the L2 convergence factor is ≤ 1/4 (+1/40).
void L2ConvergenceTest( const uint_t                                                                             minLevel,
                        const uint_t                                                                             maxLevel,
                        const System                                                                             system,
                        std::function< real_t( const uint_t level, const System& system, const bool writeVtk ) > test,
                        const bool writeVTK = false );

} // namespace n1e1
} // namespace hyteg
