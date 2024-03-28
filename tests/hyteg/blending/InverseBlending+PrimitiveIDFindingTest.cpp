/*
 * Copyright (c) 2022 Marcus Mohr.
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

// In this test we check that we can correctly reconstruct the coordinates
// of a point in the computational domain that was mapped to the physical
// domain. The catch in the test is that we assume that we do not know to
// with macro primitive the point belongs on the computational domain. This
// is e.g. the situation in the evaluate() methods of our FE function classes.

#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/math/Random.h>

#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/BlendingHelpers.hpp"
#include "hyteg/geometry/GeometryHelpers.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
// #include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

std::vector< Point3D > genSamplePointsForSimplex( uint_t dim, uint_t numSamples )
{
   walberla::math::seedRandomGenerator( 12345678 );
   std::vector< Point3D > samples;
   samples.reserve( numSamples );

   switch ( dim )
   {
   case 2: {
      while ( samples.size() < numSamples )
      {
         real_t xCoord = walberla::math::realRandom( real_c( 0 ), real_c( 1 ) );
         real_t yCoord = walberla::math::realRandom( real_c( 0 ), real_c( 1 ) );
         if ( xCoord + yCoord <= real_c( 1 ) )
         {
            Point3D sample{ { xCoord, yCoord, real_c( 0 ) } };
            samples.push_back( sample );
         }
      }
   }
   break;

   case 3: {
      while ( samples.size() < numSamples )
      {
         real_t xCoord = walberla::math::realRandom( real_c( 0 ), real_c( 1 ) );
         real_t yCoord = walberla::math::realRandom( real_c( 0 ), real_c( 1 ) );
         real_t zCoord = walberla::math::realRandom( real_c( 0 ), real_c( 1 ) );
         if ( xCoord + yCoord + zCoord <= real_c( 1 ) )
         {
            Point3D sample{ { xCoord, yCoord, zCoord } };
            samples.push_back( sample );
         }
      }
   }
   break;

   default:
      WALBERLA_ABORT( "dim = " << dim << "? You must be kidding me!" );
   }

   return samples;
}

std::vector< Point3D > genSamplePointsForFace( const Face& face, uint_t numSamples )
{
   // start by sampling reference element
   std::vector< Point3D > samples = genSamplePointsForSimplex( 2, numSamples );

   // map points to actual element in computational domain
   Matrix2r mat;
   mat( 0, 0 ) = face.getCoordinates()[1][0] - face.getCoordinates()[0][0];
   mat( 0, 1 ) = face.getCoordinates()[2][0] - face.getCoordinates()[0][0];
   mat( 1, 0 ) = face.getCoordinates()[1][1] - face.getCoordinates()[0][1];
   mat( 1, 1 ) = face.getCoordinates()[2][1] - face.getCoordinates()[0][1];
   Point2D     shift( face.getCoordinates()[0][0], face.getCoordinates()[0][1] );
   AffineMap2D affineMap( mat, shift );

   for ( auto& sample : samples )
   {
      Point3D mapped;
      affineMap.evalF( sample, mapped );
      sample = mapped;
   }

   return samples;
}

std::vector< Point3D > genSamplePointsForCell( const Cell& cell, uint_t numSamples )
{
   // start by sampling reference element
   std::vector< Point3D > samples = genSamplePointsForSimplex( 3, numSamples );

   // map points to actual element in computational domain
   Matrix3r mat;
   mat( 0, 0 ) = cell.getCoordinates()[1][0] - cell.getCoordinates()[0][0];
   mat( 0, 1 ) = cell.getCoordinates()[2][0] - cell.getCoordinates()[0][0];
   mat( 0, 2 ) = cell.getCoordinates()[3][0] - cell.getCoordinates()[0][0];

   mat( 1, 0 ) = cell.getCoordinates()[1][1] - cell.getCoordinates()[0][1];
   mat( 1, 1 ) = cell.getCoordinates()[2][1] - cell.getCoordinates()[0][1];
   mat( 1, 2 ) = cell.getCoordinates()[3][1] - cell.getCoordinates()[0][1];

   mat( 2, 0 ) = cell.getCoordinates()[1][2] - cell.getCoordinates()[0][2];
   mat( 2, 1 ) = cell.getCoordinates()[2][2] - cell.getCoordinates()[0][2];
   mat( 2, 2 ) = cell.getCoordinates()[3][2] - cell.getCoordinates()[0][2];
   AffineMap3D affineMap( mat, cell.getCoordinates()[0] );

   for ( auto& sample : samples )
   {
      Point3D mapped;
      affineMap.evalF( sample, mapped );
      sample = mapped;
   }

   return samples;
}

void testWith2DMesh( std::shared_ptr< PrimitiveStorage > storage )
{
   if ( storage->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Use testWith3DMesh() for meshes with cells!" );
   }

   // loop over all macro-faces
   for ( auto& it : storage->getFaces() )
   {
      Face& face = *it.second;

      // generate sample points on face in computational domain
      // uint_t                 numSamples      = 200;
      uint_t                 numSamples      = 20;
      std::vector< Point3D > ptsOnCompDomain = genSamplePointsForFace( face, numSamples );
      std::vector< Point3D > ptsOnPhysDomain;
      ptsOnPhysDomain.reserve( numSamples );

      // map coordinates from computational to physical domain
      for ( auto& computationalCoords : ptsOnCompDomain )
      {
         Point3D physicalCoords;
         face.getGeometryMap()->evalF( computationalCoords, physicalCoords );
         ptsOnPhysDomain.push_back( physicalCoords );
      }

      // now the actual testing
      real_t tolerance = real_c( -1 );
      if constexpr ( std::is_same_v< real_t, double > )
      {
         tolerance = real_c( 5e-15 );
      }
      else
      {
         tolerance = real_c( 5e-7 );
      }

      for ( uint_t idx = 0; idx < numSamples; ++idx )
      {
         // check that we can find the correct face for points on the computational domain
         auto [found1, faceID1] = findFaceIDForPointIn2D( storage, ptsOnCompDomain[idx], real_c( -1 ) );
         WALBERLA_CHECK( found1 );
         WALBERLA_CHECK_EQUAL( faceID1, face.getID() );

         // now try to map points back to computational domain without specifing face info
         auto [found2, faceID2, backMapped] =
             mapFromPhysicalToComputationalDomain2D( storage, ptsOnPhysDomain[idx], real_c( -1 ) );
         WALBERLA_CHECK( found2 );
         real_t reconstructionError{ ( backMapped - ptsOnCompDomain[idx] ).norm() };
         // WALBERLA_LOG_INFO_ON_ROOT( "faceID = " << std::scientific << face.getID() << ", idx = " << idx << ", error = " << reconstructionError );
         WALBERLA_CHECK_LESS_EQUAL( reconstructionError, tolerance );
         WALBERLA_CHECK_EQUAL( faceID2, face.getID() );
      }
   }
}

void testWith3DMesh( std::shared_ptr< PrimitiveStorage > storage )
{
   if ( !storage->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Use testWith2DMesh() for meshes without cells!" );
   }

   // loop over all macro-cells
   for ( auto& it : storage->getCells() )
   {
      Cell& cell = *it.second;

      // generate sample points on face in computational domain
      uint_t                 numSamples      = 10;
      std::vector< Point3D > ptsOnCompDomain = genSamplePointsForCell( cell, numSamples );
      std::vector< Point3D > ptsOnPhysDomain;
      ptsOnPhysDomain.reserve( numSamples );

      // map coordinates from computational to physical domain
      for ( auto& computationalCoords : ptsOnCompDomain )
      {
         Point3D physicalCoords;
         cell.getGeometryMap()->evalF( computationalCoords, physicalCoords );
         ptsOnPhysDomain.push_back( physicalCoords );
      }

      // now the actual testing
      real_t tolerance = real_c( -1 );
      if constexpr ( std::is_same_v< real_t, double > )
      {
         tolerance = real_c( 5e-15 );
      }
      else
      {
         tolerance = real_c( 7e-7 );
      }

      for ( uint_t idx = 0; idx < numSamples; ++idx )
      {
         // check that we can find the correct face for points on the computational domain
         auto [found1, cellID1] = findCellIDForPointIn3D( storage, ptsOnCompDomain[idx], real_c( -1 ) );
         WALBERLA_CHECK( found1 );
         WALBERLA_CHECK_EQUAL( cellID1, cell.getID() );

         // now try to map points back to computational domain without specifing face info
         auto [found2, cellID2, backMapped] =
             mapFromPhysicalToComputationalDomain3D( storage, ptsOnPhysDomain[idx], real_c( -1 ) );
         WALBERLA_CHECK( found2 );
         real_t reconstructionError{ ( backMapped - ptsOnCompDomain[idx] ).norm() };
         // WALBERLA_LOG_INFO_ON_ROOT( " idx = " << idx << std::scientific << ", error = " << reconstructionError );
         WALBERLA_CHECK_EQUAL( cellID2, cell.getID() );
         WALBERLA_CHECK_LESS_EQUAL( reconstructionError, tolerance );
         // real_t reconstructionErrorComponent0{ std::abs( backMapped[0] - ptsOnCompDomain[idx][0] ) };
         // real_t reconstructionErrorComponent1{ std::abs( backMapped[1] - ptsOnCompDomain[idx][1] ) };
         // real_t reconstructionErrorComponent2{ std::abs( backMapped[2] - ptsOnCompDomain[idx][2] ) };
         // WALBERLA_CHECK_LESS_EQUAL( reconstructionErrorComponent0, tolerance );
         // WALBERLA_CHECK_LESS_EQUAL( reconstructionErrorComponent1, tolerance );
         // WALBERLA_CHECK_LESS_EQUAL( reconstructionErrorComponent2, tolerance );
      }
   }
}

void runTests2D()
{
   // ---------
   //  Annulus
   // ---------
   real_t                rmin     = real_c( 1 );
   real_t                rmax     = real_c( 2 );
   MeshInfo              meshInfo = MeshInfo::meshAnnulus( rmin, rmax, real_c( 0 ), real_c( 2.0 * pi ), MeshInfo::CROSS, 8, 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   WALBERLA_LOG_INFO_ON_ROOT( "Running 2D test with Annulus" );
   testWith2DMesh( storage );

   // ----------------------
   //  Flow Around Cylinder
   // ----------------------
   MeshInfo              meshInfo2 = MeshInfo::fromGmshFile( "../../data/meshes/flow_around_cylinder.msh" );
   SetupPrimitiveStorage setupStorage2( meshInfo2, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   storage = std::make_shared< PrimitiveStorage >( setupStorage2 );
   WALBERLA_LOG_INFO_ON_ROOT( "Running 2D test with Flow-Around-Cylinder mesh" );
   testWith2DMesh( storage );
}

void runTests3D()
{
   // -----------------
   //  Spherical Shell
   // -----------------
   real_t                rmin     = real_c( 1 );
   real_t                rmax     = real_c( 2 );
   MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 3, 2, rmin, rmax );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   IcosahedralShellMap::setMap( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   WALBERLA_LOG_INFO_ON_ROOT( "Running 3D test with Icosahedral Shell Mesh" );
   testWith3DMesh( storage );

   // --------------------
   //  Cube with 24 Cells
   // --------------------
   MeshInfo              meshInfo2 = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_24el.msh" );
   SetupPrimitiveStorage setupStorage2( meshInfo2, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   storage = std::make_shared< PrimitiveStorage >( setupStorage2 );
   WALBERLA_LOG_INFO_ON_ROOT( "Running 3D test with cube_24el mesh" );
   testWith3DMesh( storage );
}

int main( int argc, char* argv[] )
{
   // basic setup
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // run test(s)
   runTests2D();
   runTests3D();

   return EXIT_SUCCESS;
}
