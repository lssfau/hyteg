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

#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/elementwiseoperators/N1E1ElementwiseOperator.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

void testLevel0( std::function< Eigen::Vector3r( const Point3D& ) > func, const std::array< real_t, 6 > correct )
{
   const size_t level = 0;

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   n1e1::N1E1VectorFunction< real_t > tmp( "tmp", storage, level, level );
   n1e1::N1E1VectorFunction< real_t > f( "f", storage, level, level );
   n1e1::N1E1ElementwiseMassOperator  massOp( storage, level, level );

   tmp.interpolate( func, level );
   massOp.apply( tmp, f, level, DoFType::All );

   auto cellId   = storage->getCellIDs()[0];
   auto cellData = storage->getCell( cellId )->getData( f.getDoFs()->getCellDataID() )->getPointer( level );
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yzIndex( 0, 0, 0, 0 )], correct[0], "YZ Edge (0)" )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xzIndex( 0, 0, 0, 0 )], correct[1], "XZ Edge (1)" )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xyIndex( 0, 0, 0, 0 )], correct[2], "XY Edge (2)" )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::zIndex( 0, 0, 0, 0 )], correct[3], "Z  Edge (3)" )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yIndex( 0, 0, 0, 0 )], correct[4], "Y  Edge (4)" )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xIndex( 0, 0, 0, 0 )], correct[5], "X  Edge (5)" )
}

void testLevel1()
{
   const size_t level = 1;

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   n1e1::N1E1VectorFunction< real_t > tmp( "tmp", storage, level, level );
   n1e1::N1E1VectorFunction< real_t > f( "f", storage, level, level );
   n1e1::N1E1ElementwiseMassOperator  massOp( storage, level, level );

   const Eigen::Vector3r                                    a    = { 1, 2, 3 };
   const Eigen::Vector3r                                    b    = { 4, 5, 6 };
   const std::function< Eigen::Vector3r( const Point3D& ) > func = [&]( const Point3D& x ) {
      return Eigen::Vector3r{ a + b.cross( toEigen( x ) ) };
   };

   tmp.interpolate( func, level );
   massOp.apply( tmp, f, level, DoFType::All );

   auto cellId   = storage->getCellIDs()[0];
   auto cellData = storage->getCell( cellId )->getData( f.getDoFs()->getCellDataID() )->getPointer( level );

   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xIndex( level, 0, 0, 0 )], 0.0729166666666666 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xIndex( level, 1, 0, 0 )], 0.0781250000000006 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xIndex( level, 0, 1, 0 )], 0.140625000000001 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xIndex( level, 0, 0, 1 )], 0.270833333333333 )

   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yIndex( level, 0, 0, 0 )], 0.0833333333333332 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yIndex( level, 1, 0, 0 )], 0.229166666666665 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yIndex( level, 0, 1, 0 )], 0.0729166666666644 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yIndex( level, 0, 0, 1 )], 0.135416666666666 )

   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::zIndex( level, 0, 0, 0 )], 0.0937499999999996 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::zIndex( level, 1, 0, 0 )], 0.260416666666670 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::zIndex( level, 0, 1, 0 )], 0.317708333333333 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::zIndex( level, 0, 0, 1 )], 0.0989583333333331 )

   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xyIndex( level, 0, 0, 0 )], 0.171875000000000 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xyIndex( level, 1, 0, 0 )], 0.0572916666666690 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xyIndex( level, 0, 1, 0 )], 0.0572916666666674 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xyIndex( level, 0, 0, 1 )], 0.0885416666666669 )

   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xzIndex( level, 0, 0, 0 )], 0.0104166666666666 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xzIndex( level, 1, 0, 0 )], -0.0156250000000032 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xzIndex( level, 0, 1, 0 )], 0.104166666666666 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xzIndex( level, 0, 0, 1 )], -0.0156250000000005 )

   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yzIndex( level, 0, 0, 0 )], -0.0260416666666667 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yzIndex( level, 1, 0, 0 )], -0.130208333333331 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yzIndex( level, 0, 1, 0 )], 0.0364583333333340 )
   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::yzIndex( level, 0, 0, 1 )], 0.0364583333333332 )

   WALBERLA_CHECK_FLOAT_EQUAL( cellData[edgedof::macrocell::xyzIndex( level, 0, 0, 0 )], 0.145833333333334 )
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   testLevel0( []( const Point3D& ) { return Eigen::Vector3r{ 0.0, 0.0, 0.0 }; }, { 0, 0, 0, 0, 0, 0 } );
   testLevel0(
       []( const Point3D& ) {
          return Eigen::Vector3r{ 1.0, 0.0, 0.0 };
       },
       { 0.0, -1.0 / 24.0, -1.0 / 24.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0 } );
   testLevel0(
       []( const Point3D& ) {
          return Eigen::Vector3r{ 0.0, 1.0, 0.0 };
       },
       { -1.0 / 24.0, 0.0, 1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0 } );
   testLevel0(
       []( const Point3D& ) {
          return Eigen::Vector3r{ 0.0, 0.0, 1.0 };
       },
       { 1.0 / 24.0, 1.0 / 24.0, 0.0, 1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0 } );
   testLevel0(
       []( const Point3D& p ) {
          const auto y = p[1];
          const auto z = p[2];
          return Eigen::Vector3r{ 0.0, -z, y };
       },
       { 1.0 / 30.0, 1.0 / 120.0, -1.0 / 120.0, 0.0, 0.0, 0.0 } );
   testLevel0(
       []( const Point3D& p ) {
          const auto x = p[0];
          const auto z = p[2];
          return Eigen::Vector3r{ z, 0.0, -x };
       },
       { -1.0 / 120.0, -1.0 / 30.0, -1.0 / 120.0, 0.0, 0.0, 0.0 } );
   testLevel0(
       []( const Point3D& p ) {
          const auto x = p[0];
          const auto y = p[1];
          return Eigen::Vector3r{ -y, x, 0.0 };
       },
       { -1.0 / 120.0, 1.0 / 120.0, 1.0 / 30.0, 0.0, 0.0, 0.0 } );

   testLevel1();

   return EXIT_SUCCESS;
}
