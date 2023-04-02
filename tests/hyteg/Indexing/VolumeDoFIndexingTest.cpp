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

#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"

#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/mpi/all.h"

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"
#include "hyteg/indexing/Common.hpp"

namespace hyteg {

using celldof::CellType;
using facedof::FaceType;
using indexing::Index;

/// This function tests micro-refinement index relationships.
/// Given a micro-volume index (and type) on (micro-)refinement level l, we get 4 (2D) or 8 (3D) micro-volume indices (and types)
/// on level l+1. Which ones though? This is checked here...

////////
// 2D //
////////

void testMicroRefinement2D()
{
   std::vector< Index > indexCoarse = { Index( 0, 0, 0 ), Index( 0, 0, 0 ) };

   std::vector< FaceType > typeCoarse = { FaceType::GRAY, FaceType::BLUE };

   std::vector< std::vector< Index > > indicesFine = { {
                                                           Index( 0, 0, 0 ),
                                                           Index( 1, 0, 0 ),
                                                           Index( 0, 1, 0 ),
                                                           Index( 0, 0, 0 ),
                                                       },
                                                       {
                                                           Index( 1, 0, 0 ),
                                                           Index( 0, 1, 0 ),
                                                           Index( 1, 1, 0 ),
                                                           Index( 1, 1, 0 ),

                                                       } };

   std::vector< std::vector< FaceType > > typesFine = { { FaceType::GRAY, FaceType::GRAY, FaceType::GRAY, FaceType::BLUE },
                                                        { FaceType::BLUE, FaceType::BLUE, FaceType::BLUE, FaceType::GRAY } };

   for ( uint_t i = 0; i < indexCoarse.size(); i++ )
   {
      std::vector< Index >    indexFineResult;
      std::vector< FaceType > typeFineResult;

      volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
          indexCoarse[i], typeCoarse[i], indexFineResult, typeFineResult );

      // zip all those guys to make them sortable
      std::vector< std::pair< Index, FaceType > > fineResult( 4 );
      std::vector< std::pair< Index, FaceType > > fineCheck( 4 );
      for ( size_t ii = 0; ii < 4; ii++ )
      {
         fineResult[ii].first  = indexFineResult[ii];
         fineResult[ii].second = typeFineResult[ii];

         fineCheck[ii].first  = indicesFine[i][ii];
         fineCheck[ii].second = typesFine[i][ii];
      }

      std::sort( fineResult.begin(), fineResult.end() );
      std::sort( fineCheck.begin(), fineCheck.end() );

      WALBERLA_LOG_INFO_ON_ROOT( "Result:" )
      for ( auto it : fineResult )
      {
         WALBERLA_LOG_INFO_ON_ROOT( it.first << ", " << facedof::FaceTypeToStr.at( it.second ) );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Check:" )
      for ( auto it : fineCheck )
      {
         WALBERLA_LOG_INFO_ON_ROOT( it.first << ", " << facedof::FaceTypeToStr.at( it.second ) );
      }

      WALBERLA_CHECK_EQUAL( fineResult, fineCheck );
   }
}

////////
// 3D //
////////

void testMicroRefinement3D()
{
   std::vector< Index > indexCoarse = { Index( 0, 0, 0 ) };

   std::vector< CellType > typeCoarse = { CellType::BLUE_UP };

   std::vector< std::vector< Index > > indicesFine = { { Index( 1, 0, 0 ),
                                                         Index( 0, 1, 0 ),
                                                         Index( 1, 1, 0 ),
                                                         Index( 1, 0, 1 ),
                                                         Index( 1, 0, 0 ),
                                                         Index( 1, 1, 0 ),
                                                         Index( 1, 0, 0 ),
                                                         Index( 1, 1, 0 ) } };

   std::vector< std::vector< CellType > > typesFine = { { celldof::CellType::BLUE_UP,
                                                          celldof::CellType::BLUE_UP,
                                                          celldof::CellType::BLUE_UP,
                                                          celldof::CellType::BLUE_UP,
                                                          celldof::CellType::GREEN_DOWN,
                                                          celldof::CellType::GREEN_UP,
                                                          celldof::CellType::WHITE_DOWN,
                                                          celldof::CellType::WHITE_UP } };

   for ( uint_t i = 0; i < indexCoarse.size(); i++ )
   {
      std::vector< Index >    indexFineResult;
      std::vector< CellType > typeFineResult;

      volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
          indexCoarse[i], typeCoarse[i], indexFineResult, typeFineResult );

      // zip all those guys to make them sortable
      std::vector< std::pair< Index, CellType > > fineResult( 8 );
      std::vector< std::pair< Index, CellType > > fineCheck( 8 );
      for ( size_t ii = 0; ii < 8; ii++ )
      {
         fineResult[ii].first  = indexFineResult[ii];
         fineResult[ii].second = typeFineResult[ii];

         fineCheck[ii].first  = indicesFine[i][ii];
         fineCheck[ii].second = typesFine[i][ii];
      }

      std::sort( fineResult.begin(), fineResult.end() );
      std::sort( fineCheck.begin(), fineCheck.end() );

      WALBERLA_LOG_INFO_ON_ROOT( "Result:" )
      for ( auto it : fineResult )
      {
         WALBERLA_LOG_INFO_ON_ROOT( it.first << ", " << celldof::CellTypeToStr.at( it.second ) );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "Check:" )
      for ( auto it : fineCheck )
      {
         WALBERLA_LOG_INFO_ON_ROOT( it.first << ", " << celldof::CellTypeToStr.at( it.second ) );
      }

      WALBERLA_CHECK_EQUAL( fineResult, fineCheck );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testMicroRefinement2D();
   hyteg::testMicroRefinement3D();
}
