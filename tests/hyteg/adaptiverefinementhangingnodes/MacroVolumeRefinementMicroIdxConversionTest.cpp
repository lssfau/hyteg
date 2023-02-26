/*
* Copyright (c) 2017-2022 Nils Kohl.
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
#include "core/logging/Logging.h"

#include "hyteg/Algorithms.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"

using hyteg::facedof::FaceType;
using hyteg::indexing::Index;
using walberla::int_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

/// This function tests the conversion of micro volume indices when refining volume macros.
/// There are two functions that perform these conversions - one from coarse to fine and one vice versa.
void test2D()
{
   MeshInfo              meshInfo = MeshInfo::singleTriangle( Point2D(  0, 0  ), Point2D(  1, 0  ), Point2D(  0, 1  ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   std::vector< PrimitiveID > refine = storage->getFaceIDs();
   std::vector< PrimitiveID > refineResult, coarsen, coarsenResult;

   WALBERLA_CHECK_EQUAL( refine.size(), 1 );

   storage->refinementAndCoarseningHanging( refine, coarsen, refineResult, coarsenResult );

   WALBERLA_CHECK_EQUAL( refine, refineResult );

   uint_t     level       = 3;
   const auto widthCoarse = idx_t( levelinfo::num_microedges_per_edge( level + 1 ) );

   // Test conversion from coarse idx to fine idx.
   {
      const Face& coarseFace = *storage->getFace( refine.back() );

      std::vector< Index > coarseGridIdx = { { Index( 0, 0, 0 ),
                                               Index( widthCoarse / 2, 0, 0 ),
                                               Index( 0, widthCoarse / 2, 0 ),
                                               Index( widthCoarse / 2 - 1, 0, 0 ) } };

      std::vector< FaceType > coarseGridFaceType = { { FaceType::GRAY, FaceType::GRAY, FaceType::GRAY, FaceType::BLUE } };

      std::vector< Index > fineGridIdx = { { Index( 0, 0, 0 ), Index( 0, 0, 0 ), Index( 0, 0, 0 ), Index( 0, 0, 0 ) } };

      std::vector< FaceType > fineGridFaceType = { { FaceType::GRAY, FaceType::GRAY, FaceType::GRAY, FaceType::GRAY } };

      std::vector< uint_t > fineGridMacroIdx = { { 0, 1, 2, 3 } };

      for ( uint_t i = 0; i < coarseGridIdx.size(); i++ )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "checking " << i );

         PrimitiveID finePID;
         Index       fineIdx;
         FaceType    fineFaceType;

         volumedofspace::indexing::getVolumeIdxOnRefinedMacro(
             coarseFace, level, coarseGridIdx[i], coarseGridFaceType[i], finePID, fineIdx, fineFaceType );

         WALBERLA_CHECK_EQUAL(
             fineIdx, fineGridIdx[i], "Coarse grid idx " << coarseGridIdx[i] << " gives fine grid idx: " << fineIdx );
         WALBERLA_CHECK_EQUAL( fineFaceType, fineGridFaceType[i] );
         WALBERLA_CHECK_EQUAL( finePID, coarseFace.childFaces().at( fineGridMacroIdx[i] ) );

         auto fineFace = storage->getFace( finePID );

         PrimitiveID coarsePID;
         Index       coarseIdx;
         FaceType    coarseFaceType;

         volumedofspace::indexing::getVolumeIdxOnCoarseMacro(
             *storage, *fineFace, level, fineIdx, fineFaceType, coarsePID, coarseIdx, coarseFaceType );

         WALBERLA_CHECK_EQUAL( coarsePID, coarseFace.getID() );
         WALBERLA_CHECK_EQUAL( coarseIdx, coarseGridIdx[i] );
         WALBERLA_CHECK_EQUAL( coarseFaceType, coarseGridFaceType[i] );
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   using namespace hyteg;

   test2D();
}