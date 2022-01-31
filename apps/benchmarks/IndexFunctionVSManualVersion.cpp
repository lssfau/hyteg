/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include "core/debug/TestSubsystem.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/indexing/Optimization.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

std::shared_ptr< PrimitiveStorage > globalStorage;

template < size_t Level >
inline void manualApply( real_t* oprPtr, real_t* srcPtr, real_t* dstPtr, UpdateType )
{
   size_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   size_t inner_rowsize = rowsize;

   size_t br = 1;
   size_t mr = 1 + rowsize;
   size_t tr = mr + ( rowsize - 1 );
   ///This is written for the update type Replace!
   for ( size_t i = 0; i < rowsize - 3; ++i )
   {
      for ( size_t j = 0; j < inner_rowsize - 3; ++j )
      {
         dstPtr[mr] = oprPtr[0] * srcPtr[br] + oprPtr[1] * srcPtr[br + 1] + oprPtr[2] * srcPtr[mr - 1] + oprPtr[3] * srcPtr[mr] +
                      oprPtr[4] * srcPtr[mr + 1] + oprPtr[5] * srcPtr[tr - 1] + oprPtr[6] * srcPtr[tr];

         br += 1;
         mr += 1;
         tr += 1;
      }

      br += 3;
      mr += 2;
      tr += 1;
      --inner_rowsize;
   }
}

int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
   globalStorage                               = storage;

   const size_t level = 10;

   auto                             src  = std::make_shared< hyteg::P1Function< real_t > >( "src", storage, level, level );
   auto                             dst1 = std::make_shared< hyteg::P1Function< real_t > >( "dst", storage, level, level );
   auto                             dst2 = std::make_shared< hyteg::P1Function< real_t > >( "dst", storage, level, level );
   hyteg::P1ConstantLaplaceOperator M( storage, level, level );

   std::shared_ptr< Face > face = storage->getFaces().begin().operator*().second;

   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& x ) { return x[0] * 4; };
   src->interpolate( ones, level );

   real_t* oprPtr = face->getData( M.getFaceStencilID() )->getPointer( level );
   for ( uint_t i = 0; i < 7; ++i )
   {
      oprPtr[i] = real_c( i );
   }
   real_t* srcPtr  = face->getData( src->getFaceDataID() )->getPointer( level );
   real_t* dst1Ptr = face->getData( dst1->getFaceDataID() )->getPointer( level );
   real_t* dst2Ptr = face->getData( dst2->getFaceDataID() )->getPointer( level );

   walberla::WcTimer timer;
   timer.reset();
   LIKWID_MARKER_START( "IndexApply" );
   hyteg::vertexdof::macroface::apply< real_t >(
       level, *face, M.getFaceStencilID(), src->getFaceDataID(), dst1->getFaceDataID(), Replace );
   LIKWID_MARKER_STOP( "IndexApply" );
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 20 ) << "Index Apply: " << timer.last() );

   timer.reset();
   LIKWID_MARKER_START( "ManualApply" );
   manualApply< level >( oprPtr, srcPtr, dst2Ptr, Replace );
   LIKWID_MARKER_STOP( "ManualApply" );
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( std::setw( 20 ) << "Manual Apply: " << timer.last() );

   ///check calculations
   real_t sum1 = 0, sum2 = 0;
   for ( uint_t i = 0; i < levelinfo::num_microvertices_per_face( level ); ++i )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( dst1Ptr[i], dst2Ptr[i], "i was: " << i );
      sum1 += dst1Ptr[i];
      sum2 += dst2Ptr[i];
   }
   WALBERLA_CHECK_FLOAT_UNEQUAL( sum1, 0. );
   WALBERLA_CHECK_FLOAT_UNEQUAL( sum2, 0. );

   LIKWID_MARKER_CLOSE;
}
