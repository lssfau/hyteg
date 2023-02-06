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

#include "hyteg/gridtransferoperators/P1toN1E1Gradient.hpp"

#include <Eigen/Sparse>

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

void test( const uint_t                                              lvl,
           const std::function< real_t( const Point3D& ) >&          f,
           const std::function< Eigen::Vector3r( const Point3D& ) >& gradF,
           const bool                                                writeVTK = false )
{
   using namespace n1e1;

   MeshInfo              meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P1Function< real_t > p1F( "p1F", storage, lvl, lvl );
   p1F.interpolate( f, lvl );

   N1E1VectorFunction< real_t > n1e1F( "n1e1F", storage, lvl, lvl );
   n1e1F.interpolate( Eigen::Vector3r{ -1.1, 1000.0, 4.2 }, lvl ); // junk

   N1E1VectorFunction< real_t > projectedAnalytical( "projected analytical", storage, lvl, lvl );
   projectedAnalytical.interpolate( gradF, lvl );
   projectedAnalytical.communicate< Edge, Face >( lvl );

   P1toN1E1Gradient( p1F, n1e1F, lvl );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output/", "P1toN1E1GradientTest", storage );
      vtk.add( p1F );
      vtk.add( n1e1F );
      vtk.add( projectedAnalytical );
      vtk.write( lvl );
   }

   for ( const auto& it : storage->getEdges() )
   {
      using edgedof::macroedge::index;

      const Edge& edge        = *it.second;
      const auto  testData    = edge.getData( n1e1F.getDoFs()->getEdgeDataID() )->getPointer( lvl );
      const auto  correctData = edge.getData( projectedAnalytical.getDoFs()->getEdgeDataID() )->getPointer( lvl );

      for ( auto idx : edgedof::macroedge::Iterator( lvl ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( testData[index( lvl, idx.x() )], correctData[index( lvl, idx.x() )] )
      }
   }

   for ( const auto& it : storage->getFaces() )
   {
      using edgedof::macroface::diagonalIndex;
      using edgedof::macroface::horizontalIndex;
      using edgedof::macroface::verticalIndex;

      const Face& face        = *it.second;
      const auto  testData    = face.getData( n1e1F.getDoFs()->getFaceDataID() )->getPointer( lvl );
      const auto  correctData = face.getData( projectedAnalytical.getDoFs()->getFaceDataID() )->getPointer( lvl );

      for ( auto idx : edgedof::macroface::Iterator( lvl ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( testData[horizontalIndex( lvl, idx.x(), idx.y() )],
                                     correctData[horizontalIndex( lvl, idx.x(), idx.y() )] )
         WALBERLA_CHECK_FLOAT_EQUAL( testData[verticalIndex( lvl, idx.x(), idx.y() )],
                                     correctData[verticalIndex( lvl, idx.x(), idx.y() )] )
         WALBERLA_CHECK_FLOAT_EQUAL( testData[diagonalIndex( lvl, idx.x(), idx.y() )],
                                     correctData[diagonalIndex( lvl, idx.x(), idx.y() )] )
      }
   }

   for ( const auto& it : storage->getCells() )
   {
      using edgedof::macrocell::xIndex;
      using edgedof::macrocell::xyIndex;
      using edgedof::macrocell::xyzIndex;
      using edgedof::macrocell::xzIndex;
      using edgedof::macrocell::yIndex;
      using edgedof::macrocell::yzIndex;
      using edgedof::macrocell::zIndex;

      const Cell& cell        = *it.second;
      const auto  testData    = cell.getData( n1e1F.getDoFs()->getCellDataID() )->getPointer( lvl );
      const auto  correctData = cell.getData( projectedAnalytical.getDoFs()->getCellDataID() )->getPointer( lvl );

      for ( auto idx : edgedof::macrocell::Iterator( lvl ) )
      {
         const idx_t x = idx.x();
         const idx_t y = idx.y();
         const idx_t z = idx.z();

         WALBERLA_CHECK_FLOAT_EQUAL( testData[xIndex( lvl, x, y, z )], correctData[xIndex( lvl, x, y, z )] )
         WALBERLA_CHECK_FLOAT_EQUAL( testData[yIndex( lvl, x, y, z )], correctData[yIndex( lvl, x, y, z )] )
         WALBERLA_CHECK_FLOAT_EQUAL( testData[zIndex( lvl, x, y, z )], correctData[zIndex( lvl, x, y, z )] )
         WALBERLA_CHECK_FLOAT_EQUAL( testData[xyIndex( lvl, x, y, z )], correctData[xyIndex( lvl, x, y, z )] )
         WALBERLA_CHECK_FLOAT_EQUAL( testData[xzIndex( lvl, x, y, z )], correctData[xzIndex( lvl, x, y, z )] )
         WALBERLA_CHECK_FLOAT_EQUAL( testData[yzIndex( lvl, x, y, z )], correctData[yzIndex( lvl, x, y, z )] )
      }

      for ( auto idx : edgedof::macrocell::IteratorXYZ( lvl ) )
      {
         const idx_t x = idx.x();
         const idx_t y = idx.y();
         const idx_t z = idx.z();

         WALBERLA_CHECK_FLOAT_EQUAL( testData[xyzIndex( lvl, x, y, z )], correctData[xyzIndex( lvl, x, y, z )] )
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::function< real_t( const Point3D& ) >          f     = []( const Point3D& ) { return 42.0; };
   std::function< Eigen::Vector3r( const Point3D& ) > gradF = []( const Point3D& ) { return Eigen::Vector3r{ 0.0, 0.0, 0.0 }; };

   std::function< real_t( const Point3D& ) >          g = []( const Point3D& p ) { return 3.0 * p[0] + 5.0 * p[1] + 7.0 * p[2]; };
   std::function< Eigen::Vector3r( const Point3D& ) > gradG = []( const Point3D& ) { return Eigen::Vector3r{ 3.0, 5.0, 7.0 }; };

   // Note that the following functions pass the test with exact, that is independent of h, checks
   // although neither function is in its respective FEM-space.
   // This is because the diagram below commutes (we hereby test that this is true):
   //
   //        grad
   //    H¹ -----→ H(curl)
   //    |         |
   // Π⁰ |         | Π¹
   //    ↓   grad  ↓
   //    P1 -----→ N1E1
   //
   // Note also that this only works with linear gradients because we use the mid-point rule to approximate Π¹.
   std::function< real_t( const Point3D& ) > h = []( const Point3D& p ) {
      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];
      return 0.5 * x * x - y * y + 1.5 * z * z + 0.5 * x * y - 4.0 * x * z + 2.5 * y * z;
   };
   std::function< Eigen::Vector3r( const Point3D& ) > gradH = []( const Point3D& p ) {
      const real_t x = p[0];
      const real_t y = p[1];
      const real_t z = p[2];
      return Eigen::Vector3r{
          real_c( x + 0.5 * y - 4.0 * z ), real_c( 0.5 * x - 2.0 * y + 2.5 * z ), real_c( -4.0 * x + 2.5 * y + 3.0 * z ) };
   };

   test( 1, f, gradF );
   test( 3, g, gradG );
   test( 3, h, gradH );

   return EXIT_SUCCESS;
}
