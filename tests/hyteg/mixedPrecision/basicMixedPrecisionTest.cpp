/*
* Copyright (c) 2017-2023 Dominik Thoennes.
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
#include "core/debug/Debug.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

void testMixedPrecision()
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t                                           minLevel   = 4;
   uint_t                                           maxLevel   = 4;
   std::function< double( const hyteg::Point3D& ) > doubleRand = []( const hyteg::Point3D& ) {
      return walberla::math::realRandom< double >();
   };
   std::function< float( const hyteg::Point3D& ) > floatRand = []( const hyteg::Point3D& ) {
      return walberla::math::realRandom< float >();
   };
   auto p1Double = P1Function< double >( "doubleP1", storage, minLevel, maxLevel );
   auto p1Float  = P1Function< float >( "doubleP1", storage, minLevel, maxLevel );

   p1Double.interpolate( doubleRand, minLevel );
   p1Float.copyFrom( p1Double, minLevel );

   auto doubleFaceData  = storage->getFace( storage->getFaceIDs()[0] )->getData( p1Double.getFaceDataID() );
   auto p1DobuleFacePtr = doubleFaceData->getPointer( minLevel );
   auto floatFaceData   = storage->getFace( storage->getFaceIDs()[0] )->getData( p1Float.getFaceDataID() );
   auto p1FloatFacePtr  = floatFaceData->getPointer( minLevel );

   for ( uint_t i = 0; i < doubleFaceData->getSize( minLevel ); ++i )
   {
      // The difference between float and double should not be bigger than 1e-7
      WALBERLA_CHECK_LESS( p1DobuleFacePtr[i] - static_cast< double >( p1FloatFacePtr[i] ), 1e-7 )
   }

   p1Float.interpolate( floatRand, minLevel );
   p1Double.copyFrom( p1Float, minLevel );

   for ( uint_t i = 0; i < doubleFaceData->getSize( minLevel ); ++i )
   {
      WALBERLA_CHECK_EQUAL( p1DobuleFacePtr[i] - static_cast< double >( p1FloatFacePtr[i] ), real_t( 0 ) )
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testMixedPrecision();

   return EXIT_SUCCESS;
}
