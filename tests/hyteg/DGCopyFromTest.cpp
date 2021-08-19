/*
 * Copyright (c) 2021 Andreas Wagner.
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
#include "core/debug/all.h"
#include "core/mpi/all.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/facedofspace/FaceDoFFunction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using namespace hyteg;

using walberla::real_t;

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::debug::enterTestMode();

   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage1 = std::make_shared< PrimitiveStorage >( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage2 = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   hyteg::FaceDoFFunction< real_t > x1( "x1", storage1, minLevel, maxLevel );
   hyteg::FaceDoFFunction< real_t > x2( "x2", storage2, minLevel, maxLevel );

   x1.interpolate(1, maxLevel, All);
   x2.copyFrom(x1, maxLevel);

   // check if copying worked:
   WALBERLA_CHECK_FLOAT_EQUAL( x2.getMinValue( maxLevel, All ), 1. );
   WALBERLA_CHECK_FLOAT_EQUAL( x2.getMaxValue( maxLevel, All ), 1. );

   // check that it is a copy and changing x1 has no impact on x2:
   x1.interpolate(0, maxLevel, All);

   WALBERLA_CHECK_FLOAT_EQUAL( x2.getMinValue( maxLevel, All ), 1. );
   WALBERLA_CHECK_FLOAT_EQUAL( x2.getMaxValue( maxLevel, All ), 1. );

   WALBERLA_CHECK_FLOAT_EQUAL( x1.getMinValue( maxLevel, All ), 0. );
   WALBERLA_CHECK_FLOAT_EQUAL( x1.getMaxValue( maxLevel, All ), 0. );

   return 0;
}
