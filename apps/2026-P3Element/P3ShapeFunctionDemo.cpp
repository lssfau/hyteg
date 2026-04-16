/*
 * Copyright (c) 2026 Marcus Mohr.
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
#include "core/Environment.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p3functionspace/P3Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;

using triplet = std::tuple< real_t, real_t, real_t >;

void printSeparator()
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------------" );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "Shape Functions" );

   // using the unit simplex for this demo
   MeshInfo              mesh2D = MeshInfo::singleTriangle( Point2D{ 0.0, 0.0 }, Point2D{ 1.0, 0.0 }, Point2D{ 0.0, 1.0 } );
   SetupPrimitiveStorage setupStorage( mesh2D, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t level = 5;

   // auxilliary function for computing barycentric coordinates of a point
   std::function< triplet( const hyteg::Point3D& ) > getBarycentricCoordinates = []( const Point3D& x ) -> triplet {
      real_t L1 = real_c( 1 ) - x[0] - x[1];
      real_t L2 = x[0];
      real_t L3 = x[1];
      return { L1, L2, L3 };
   };

   // -------------
   //  Vertex DoFs
   // -------------

   // shape function associated with vertex #0
   std::function< real_t( const hyteg::Point3D& ) > shape0 = []( const Point3D& x ) -> real_t {
      real_t L1 = real_c( 1 ) - x[0] - x[1];
      return real_c( 0.5 ) * L1 * ( real_c( 3 ) * L1 - real_c( 1 ) ) * ( real_c( 3 ) * L1 - real_c( 2 ) );
   };
   P3Function< real_t > shapeFunc0( "Shape Function 0", storage, level, level );
   shapeFunc0.interpolate( shape0, level, All );

   // shape function associated with vertex #1
   std::function< real_t( const hyteg::Point3D& ) > shape1 = []( const Point3D& x ) -> real_t {
      real_t L2 = x[0];
      return real_c( 0.5 ) * L2 * ( real_c( 3 ) * L2 - real_c( 1 ) ) * ( real_c( 3 ) * L2 - real_c( 2 ) );
   };
   P3Function< real_t > shapeFunc1( "Shape Function 1", storage, level, level );
   shapeFunc1.interpolate( shape1, level, All );

   // shape function associated with vertex #2
   std::function< real_t( const hyteg::Point3D& ) > shape2 = []( const Point3D& x ) -> real_t {
      real_t L3 = x[1];
      return real_c( 0.5 ) * L3 * ( real_c( 3 ) * L3 - real_c( 1 ) ) * ( real_c( 3 ) * L3 - real_c( 2 ) );
   };
   P3Function< real_t > shapeFunc2( "Shape Function 2", storage, level, level );
   shapeFunc2.interpolate( shape2, level, All );

   // -----------
   //  Edge DoFs
   // -----------

   // shape function associated with vertex #3 (blue DoF on edge 1 -> 2)
   std::function< real_t( const hyteg::Point3D& ) > shape3 = []( const Point3D& x ) -> real_t {
      real_t L2 = x[0];
      real_t L3 = x[1];
      return real_c( 4.5 ) * L2 * L3 * ( real_c( 3 ) * L2 - real_c( 1 ) );
   };
   P3Function< real_t > shapeFunc3( "Shape Function 3", storage, level, level );
   shapeFunc3.interpolate( shape3, level, All );

   // shape function associated with vertex #4 (blue DoF on edge 0 -> 2)
   std::function< real_t( const hyteg::Point3D& ) > shape4 = []( const Point3D& x ) -> real_t {
      real_t L1 = real_c( 1 ) - x[0] - x[1];
      real_t L3 = x[1];
      return real_c( 4.5 ) * L3 * L1 * ( real_c( 3 ) * L1 - real_c( 1 ) );
   };
   P3Function< real_t > shapeFunc4( "Shape Function 4", storage, level, level );
   shapeFunc4.interpolate( shape4, level, All );

   // shape function associated with vertex #5 (blue DoF on edge 0 -> 1)
   std::function< real_t( const hyteg::Point3D& ) > shape5 = []( const Point3D& x ) -> real_t {
      real_t L1 = real_c( 1 ) - x[0] - x[1];
      real_t L2 = x[0];
      return real_c( 4.5 ) * L1 * L2 * ( real_c( 3 ) * L1 - real_c( 1 ) );
   };
   P3Function< real_t > shapeFunc5( "Shape Function 5", storage, level, level );
   shapeFunc5.interpolate( shape5, level, All );

   // shape function associated with vertex #6 (red DoF on edge 1 -> 2)
   std::function< real_t( const hyteg::Point3D& ) > shape6 = []( const Point3D& x ) -> real_t {
      real_t L2 = x[0];
      real_t L3 = x[1];
      return real_c( 4.5 ) * L2 * L3 * ( real_c( 3 ) * L3 - real_c( 1 ) );
   };
   P3Function< real_t > shapeFunc6( "Shape Function 6", storage, level, level );
   shapeFunc6.interpolate( shape6, level, All );

   // shape function associated with vertex #7 (red DoF on edge 0 -> 2)
   std::function< real_t( const hyteg::Point3D& ) > shape7 = []( const Point3D& x ) -> real_t {
      real_t L1 = real_c( 1 ) - x[0] - x[1];
      real_t L3 = x[1];
      return real_c( 4.5 ) * L3 * L1 * ( real_c( 3 ) * L3 - real_c( 1 ) );
   };
   P3Function< real_t > shapeFunc7( "Shape Function 7", storage, level, level );
   shapeFunc7.interpolate( shape7, level, All );

   // shape function associated with vertex #8 (red DoF on edge 0 -> 1)
   std::function< real_t( const hyteg::Point3D& ) > shape8 = []( const Point3D& x ) -> real_t {
      real_t L1 = real_c( 1 ) - x[0] - x[1];
      real_t L2 = x[0];
      return real_c( 4.5 ) * L1 * L2 * ( real_c( 3 ) * L2 - real_c( 1 ) );
   };
   P3Function< real_t > shapeFunc8( "Shape Function 8", storage, level, level );
   shapeFunc8.interpolate( shape8, level, All );

   // ----------
   //  Face DoF
   // ----------
 
   // shape function associated with vertex #9 (barycenter)
   std::function< real_t( const hyteg::Point3D& ) > shape9 = []( const Point3D& x ) -> real_t {
      real_t L1 = real_c( 1 ) - x[0] - x[1];
      real_t L2 = x[0];
      real_t L3 = x[1];
      return real_c( 27 ) * L1 * L2 * L3;
   };
   P3Function< real_t > shapeFunc9( "Shape Function 9", storage, level, level );
   shapeFunc9.interpolate( shape9, level, All );

   // --------
   //  Export
   // ---------

   VTKOutput vtkOutput( "output", "P3ShapeFunctions", storage );
   vtkOutput.add( shapeFunc0 );
   vtkOutput.add( shapeFunc1 );
   vtkOutput.add( shapeFunc2 );
   vtkOutput.add( shapeFunc3 );
   vtkOutput.add( shapeFunc4 );
   vtkOutput.add( shapeFunc5 );
   vtkOutput.add( shapeFunc6 );
   vtkOutput.add( shapeFunc7 );
   vtkOutput.add( shapeFunc8 );
   vtkOutput.add( shapeFunc9 );
   vtkOutput.write( level );

   return EXIT_SUCCESS;
}
