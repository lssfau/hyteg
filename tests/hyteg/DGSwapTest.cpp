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
#include "hyteg/facedofspace_old/FaceDoFFunction.hpp"
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
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   hyteg::FaceDoFFunction_old< real_t > x( "x", storage, minLevel, maxLevel );
   hyteg::FaceDoFFunction_old< real_t > y( "y", storage, minLevel, maxLevel );

   x.interpolate( 1, maxLevel, hyteg::All );
   y.interpolate( 2, maxLevel, hyteg::All );

   x.swap( y, maxLevel, All );

   // check if swapping was successful:
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, All ), 2. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, All ), 2. );

   WALBERLA_CHECK_FLOAT_EQUAL( y.getMinValue( maxLevel, All ), 1. );
   WALBERLA_CHECK_FLOAT_EQUAL( y.getMaxValue( maxLevel, All ), 1. );

   // check if changing x modifies y, i.e. if parts were not swapped:
   x.interpolate( 3, maxLevel, All );

   y.communicate< Vertex, Edge >( maxLevel );
   y.communicate< Edge, Face >( maxLevel );

   WALBERLA_CHECK_FLOAT_EQUAL( y.getMinValue( maxLevel, All ), 1. );
   WALBERLA_CHECK_FLOAT_EQUAL( y.getMaxValue( maxLevel, All ), 1. );

   return 0;
}
