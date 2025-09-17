/*
 * Copyright (c) 2025 Ponsuganth Ilangovan
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

// #include "terraneo/utils/MeshTransferWithParticles.hpp"

#include "terraneo/utils/MeshTransferWithParticles.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg_operators/operators/mass/P1ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassAnnulusMap.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"

using namespace hyteg;

template < typename FunctionType, typename BlendingMapType, typename MassOperatorType >
void testMeshTransfer( MeshInfo meshInfoSrc, MeshInfo meshInfoDst, uint_t levelSrc, uint_t levelDst, real_t errorAllowed )
{
   SetupPrimitiveStorage setupStorage( meshInfoSrc, walberla::MPIManager::instance()->numProcesses() );
   BlendingMapType::setMap( setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1u );
   storage->keepAllPrimitivesAsNeighbors( setupStorage );

   FunctionType T( "T", storage, levelSrc, levelSrc );
   FunctionType test( "test", storage, levelSrc, levelSrc );

   std::function< real_t( const Point3D& ) > TInterp = []( const Point3D& X ) {
      real_t x = X[0];
      real_t y = X[1];
      real_t z = X[2];
      return std::sin( x ) * std::sin( x ) + y * y * std::cos( z );
   };

   T.interpolate( TInterp, levelSrc, All );
   communication::syncFunctionBetweenPrimitives( T, levelSrc );

   SetupPrimitiveStorage setupStorageNew( meshInfoDst, walberla::MPIManager::instance()->numProcesses() );
   BlendingMapType::setMap( setupStorageNew );
   auto storageNew = std::make_shared< PrimitiveStorage >( setupStorageNew, 1u );

   FunctionType TNew( "TNew", storageNew, levelDst, levelDst );
   FunctionType testNew( "testNew", storageNew, levelDst, levelDst );

   auto meshTransferUtility = MeshTransferWithParticles< FunctionType >();
   meshTransferUtility.transfer( T, TNew, levelSrc, levelDst );

   MassOperatorType massOperator( storage, levelSrc, levelSrc );
   MassOperatorType massOperatorNew( storageNew, levelDst, levelDst );

   massOperator.apply( T, test, levelSrc, All );
   massOperatorNew.apply( TNew, testNew, levelDst, All );

   real_t testSum    = test.sumGlobal( levelSrc, All );
   real_t testNewSum = testNew.sumGlobal( levelDst, All );
   real_t difference = std::abs( testSum - testNewSum );

   // WALBERLA_LOG_INFO_ON_ROOT( "testSum = " << testSum );
   // WALBERLA_LOG_INFO_ON_ROOT( "testNewSum = " << testNewSum );
   // WALBERLA_LOG_INFO_ON_ROOT( "difference = " << std::abs( testSum - testNewSum ) );

   WALBERLA_CHECK_LESS( difference, errorAllowed );
}

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo meshInfoSquareSrc = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, 4u, 4u );
   MeshInfo meshInfoSquareDst = MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, 8u, 8u );

   testMeshTransfer< P1Function< real_t >, IdentityMap, operatorgeneration::P1ElementwiseMass >(
       meshInfoSquareSrc, meshInfoSquareDst, 2u, 1u, 1e-5 );
   testMeshTransfer< P2Function< real_t >, IdentityMap, operatorgeneration::P2ElementwiseMass >(
       meshInfoSquareSrc, meshInfoSquareDst, 2u, 1u, 1e-5 );

   MeshInfo meshInfoAnnulusSrc = MeshInfo::meshAnnulus( 1.22, 2.22, MeshInfo::CRISS, 4u, 2u );
   MeshInfo meshInfoAnnulusDst = MeshInfo::meshAnnulus( 1.22, 2.22, MeshInfo::CRISS, 8u, 4u );

   testMeshTransfer< P2Function< real_t >, AnnulusMap, operatorgeneration::P2ElementwiseMassAnnulusMap >(
       meshInfoAnnulusSrc, meshInfoAnnulusDst, 3u, 2u, 1e-2 );

   MeshInfo meshInfoCuboidSrc = MeshInfo::meshCuboid( Point3D( 0.0, 0.0, 0.0 ), Point3D( 1.0, 1.0, 1.0 ), 1u, 1u, 1u );
   MeshInfo meshInfoCuboidDst = MeshInfo::meshCuboid( Point3D( 0.0, 0.0, 0.0 ), Point3D( 1.0, 1.0, 1.0 ), 2u, 2u, 2u );

   testMeshTransfer< P1Function< real_t >, IdentityMap, operatorgeneration::P1ElementwiseMass >(
       meshInfoCuboidSrc, meshInfoCuboidDst, 2u, 1u, 1e-2 );
   testMeshTransfer< P2Function< real_t >, IdentityMap, operatorgeneration::P2ElementwiseMass >(
       meshInfoCuboidSrc, meshInfoCuboidDst, 2u, 1u, 1e-2 );

   MeshInfo meshInfoSphSrc = MeshInfo::meshSphericalShell( 3u, 2u, 1.22, 2.22 );
   MeshInfo meshInfoSphDst = MeshInfo::meshSphericalShell( 3u, 3u, 1.22, 2.22 );

   testMeshTransfer< P2Function< real_t >, IcosahedralShellMap, operatorgeneration::P2ElementwiseMassIcosahedralShellMap >(
       meshInfoSphSrc, meshInfoSphDst, 2u, 1u, 1e-2 );

   return EXIT_SUCCESS;
}