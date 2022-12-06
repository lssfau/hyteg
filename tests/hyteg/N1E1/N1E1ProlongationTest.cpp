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
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/eigen/typeAliases.hpp"
#include "hyteg/gridtransferoperators/N1E1toN1E1Prolongation.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

void test3D()
{
   MeshInfo              meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const size_t minLevel    = 3;
   const size_t maxLevel    = 4;
   const size_t coarseLevel = maxLevel - 1;
   const size_t fineLevel   = maxLevel;

   walberla::math::seedRandomGenerator( 42 );

   const uint_t numRandomEvaluations = 1000;

   // most general function in N1E1 space
   const Eigen::Vector3r                                    a          = { 1, 2, 3 };
   const Eigen::Vector3r                                    b          = { 4, 5, 6 };
   const std::function< Eigen::Vector3r( const Point3D& ) > testFuncAB = [&]( const Point3D& x ) {
      return ( a + b.cross( x.vector_ ) ).eval();
   };
   const std::function< Eigen::Vector3r( const Point3D& ) > testFuncBA = [&]( const Point3D& x ) {
      return ( b + a.cross( x.vector_ ) ).eval();
   };

   n1e1::N1E1VectorFunction< real_t > f( "f", storage, minLevel, maxLevel );
   f.interpolate( testFuncAB, coarseLevel );
   f.interpolate( testFuncBA, fineLevel );

   n1e1::N1E1VectorFunction< real_t > g( "g", storage, minLevel, maxLevel );
   g.interpolate( testFuncAB, coarseLevel );
   g.interpolate( testFuncBA, fineLevel );

   n1e1::N1E1VectorFunction< real_t > ba( "ba", storage, minLevel, maxLevel );
   ba.interpolate( testFuncBA, coarseLevel );

   n1e1::N1E1toN1E1Prolongation prolongation;
   prolongation.prolongate( f, coarseLevel, All );
   prolongation.prolongateAndAdd( g, coarseLevel, All );

   // test on cells
   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      Point3D coordinates;
      coordinates[0] = walberla::math::realRandom( 0.0, 1.0 );
      coordinates[1] = walberla::math::realRandom( 0.0, 1.0 );
      coordinates[2] = walberla::math::realRandom( 0.0, 1.0 );

      Eigen::Vector3r evalCoarse, evalFine, evalBA;
      auto            successCoarse = f.evaluate( coordinates, coarseLevel, evalCoarse );
      auto            successFine   = f.evaluate( coordinates, fineLevel, evalFine );
      WALBERLA_CHECK( successCoarse );
      WALBERLA_CHECK( successFine );
      WALBERLA_CHECK_FLOAT_EQUAL( evalCoarse[0], evalFine[0], "Test3D: wrong X-coordinate at " << coordinates << "." );
      WALBERLA_CHECK_FLOAT_EQUAL( evalCoarse[1], evalFine[1], "Test3D: wrong Y-coordinate at " << coordinates << "." );
      WALBERLA_CHECK_FLOAT_EQUAL( evalCoarse[2], evalFine[2], "Test3D: wrong Z-coordinate at " << coordinates << "." );

      successFine          = g.evaluate( coordinates, fineLevel, evalFine );
      const bool successBA = ba.evaluate( coordinates, coarseLevel, evalBA );
      WALBERLA_CHECK( successFine );
      WALBERLA_CHECK( successBA );
      WALBERLA_CHECK_FLOAT_EQUAL(
          evalCoarse[0] + evalBA[0], evalFine[0], "Additive test3D: wrong X-coordinate at " << coordinates << "." );
      WALBERLA_CHECK_FLOAT_EQUAL(
          evalCoarse[1] + evalBA[1], evalFine[1], "Additive test3D: wrong Y-coordinate at " << coordinates << "." );
      WALBERLA_CHECK_FLOAT_EQUAL(
          evalCoarse[2] + evalBA[2], evalFine[2], "Additive test3D: wrong Z-coordinate at " << coordinates << "." );
   }

   // test on edges and faces
   n1e1::N1E1VectorFunction< real_t > tmpF( "tmpF", storage, minLevel, maxLevel );
   n1e1::N1E1VectorFunction< real_t > tmpG( "tmpG", storage, minLevel, maxLevel );
   tmpF.copyFrom( f, fineLevel );
   tmpG.copyFrom( g, fineLevel );

   tmpF.communicate< Edge, Face >( fineLevel );
   tmpF.communicate< Face, Cell >( fineLevel );

   tmpG.communicate< Edge, Face >( fineLevel );
   tmpG.communicate< Face, Cell >( fineLevel );

   for ( auto& it : storage->getCells() )
   {
      Cell&      cell     = *it.second;
      const auto fData    = cell.getData( f.getDoFs()->getCellDataID() )->getPointer( fineLevel );
      const auto tmpFData = cell.getData( tmpF.getDoFs()->getCellDataID() )->getPointer( fineLevel );

      const auto gData    = cell.getData( g.getDoFs()->getCellDataID() )->getPointer( fineLevel );
      const auto tmpGData = cell.getData( tmpG.getDoFs()->getCellDataID() )->getPointer( fineLevel );

      for ( auto idx : edgedof::macrocell::Iterator( fineLevel ) )
      {
         const idx_t x = idx.x();
         const idx_t y = idx.y();
         const idx_t z = idx.z();

         for ( const auto edgeType : edgedof::allEdgeDoFOrientationsWithoutXYZ )
         {
            const real_t fVal    = fData[edgedof::macrocell::index( fineLevel, x, y, z, edgeType )];
            const real_t tmpFVal = tmpFData[edgedof::macrocell::index( fineLevel, x, y, z, edgeType )];
            WALBERLA_CHECK_FLOAT_EQUAL( fVal,
                                        tmpFVal,
                                        "edge type = " << edgeType << ", idx = " << idx << ", cell coordinates = ["
                                                       << cell.getCoordinates()[0] << ", " << cell.getCoordinates()[1] << ", "
                                                       << cell.getCoordinates()[2] << ", " << cell.getCoordinates()[3] << "]" );

            const real_t gVal    = gData[edgedof::macrocell::index( fineLevel, x, y, z, edgeType )];
            const real_t tmpGVal = tmpGData[edgedof::macrocell::index( fineLevel, x, y, z, edgeType )];
            WALBERLA_CHECK_FLOAT_EQUAL( gVal,
                                        tmpGVal,
                                        "edge type = " << edgeType << ", idx = " << idx << ", cell coordinates = ["
                                                       << cell.getCoordinates()[0] << ", " << cell.getCoordinates()[1] << ", "
                                                       << cell.getCoordinates()[2] << ", " << cell.getCoordinates()[3] << "]" );
         }
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test3D();

   return EXIT_SUCCESS;
}
