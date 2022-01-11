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

   x.interpolate( -1, maxLevel, DirichletBoundary );
   x.interpolate( 2, maxLevel, Inner );

   WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, DirichletBoundary ), -1. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, DirichletBoundary ), -1. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxMagnitude( maxLevel, DirichletBoundary ), 1. );

   WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, Inner ), 2. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, Inner ), 2. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxMagnitude( maxLevel, Inner ), 2. );

   WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, All ), -1. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, All ), 2. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxMagnitude( maxLevel, All ), 2. );

   const double h = 1. / levelinfo::num_microedges_per_edge( maxLevel );

   x.interpolate( []( const auto& p ) { return p[0] - 0.5; }, maxLevel, All );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMinValue( maxLevel, All ), -0.5 + h / 3. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxValue( maxLevel, All ), 0.5 - h / 3. );
   WALBERLA_CHECK_FLOAT_EQUAL( x.getMaxMagnitude( maxLevel, All ), 0.5 - h / 3. );

   return 0;
}
