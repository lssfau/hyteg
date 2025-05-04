/*
 * Copyright (c) 2025 Andreas Burkhart.
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

#include <cmath>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/indexing/ConsistentEnumeration.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/sparseassembly/FileWritingVector.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   uint_t minLevel_ = 0;
   uint_t maxLevel_ = 3;
   uint_t n_        = 2;

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISSCROSS, n_, n_ );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // set load balancing
   loadbalancing::roundRobinVolume( setupStorage );

   // set boundary flags
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // create storage
   auto storage_ = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

   BoundaryCondition bc = BoundaryCondition::create0123BC();

   auto T1_        = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "T1", storage_, minLevel_, maxLevel_, bc );
   auto T2_        = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "T2", storage_, minLevel_, maxLevel_, bc );
   auto T3_        = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "T3", storage_, minLevel_, maxLevel_, bc );
   auto enumerator = std::make_shared< P2P1TaylorHoodFunction< idx_t > >( "enumerator", storage_, minLevel_, maxLevel_, bc );

   std::function< real_t( const hyteg::Point3D& ) > k = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return x[1] * x[0] - x[0] - 3;
   };

   T1_->interpolate( k, maxLevel_, All );

   WALBERLA_LOG_INFO_ON_ROOT( "Start consistent enum..." );
   enumerateConsistently( *enumerator, maxLevel_ );
   WALBERLA_LOG_INFO_ON_ROOT( "finish." );

   auto vec  = std::make_shared< FileWritingVector< P2P1TaylorHoodFunction< real_t > > >( storage_, maxLevel_, *enumerator, T1_ );
   auto vec2 = std::make_shared< FileWritingVector< P2P1TaylorHoodFunction< real_t > > >( storage_, maxLevel_, *enumerator, T2_ );

   WALBERLA_LOG_INFO_ON_ROOT( "Start writing to file..." );
   vec->writeToFile( "FileWritingVectorAndConsistentEnumerationTest.bin" );
   WALBERLA_LOG_INFO_ON_ROOT( "finish." );

   WALBERLA_LOG_INFO_ON_ROOT( "Start loading from file..." );
   vec2->readFromFile( "FileWritingVectorAndConsistentEnumerationTest.bin" );
   WALBERLA_LOG_INFO_ON_ROOT( "finish." );

   T2_->fromVector( *enumerator, vec2, maxLevel_, All );

   T3_->assign( { real_c( 1 ), real_c( -1 ) }, { *T1_, *T2_ }, maxLevel_, All );

   real_t     DiffSum    = T3_->dotGlobal( *T3_, maxLevel_, All );
   const bool DiffSumIsZero = std::fpclassify( DiffSum ) == FP_ZERO;
   WALBERLA_LOG_INFO_ON_ROOT( "Load/Reload DiffSum = " << DiffSum );
   WALBERLA_CHECK( DiffSumIsZero );
}