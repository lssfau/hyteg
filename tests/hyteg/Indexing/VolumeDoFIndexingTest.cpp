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
   const uint_t levelCoarse = 3;

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

      auto indexFineCheck = indicesFine[i];
      auto typeFineCheck  = typesFine[i];

      auto permutationCheck  = algorithms::sortPermutation( indexFineCheck, std::less< Index >() );
      auto permutationResult = algorithms::sortPermutation( indexFineResult, std::less< Index >() );

      indexFineCheck = algorithms::applyPermutation( indexFineCheck, permutationCheck );
      typeFineCheck  = algorithms::applyPermutation( typeFineCheck, permutationCheck );

      indexFineResult = algorithms::applyPermutation( indexFineResult, permutationResult );
      typeFineResult  = algorithms::applyPermutation( typeFineResult, permutationResult );

      WALBERLA_LOG_INFO_ON_ROOT( "checking coarse index: " << indexCoarse[i]
                                                           << ", type: " << facedof::FaceTypeToStr.at( typeCoarse[i] ) );

      WALBERLA_CHECK_EQUAL( indexFineCheck, indexFineResult );
      WALBERLA_CHECK_EQUAL( typeFineCheck, typeFineResult );
   }
}

////////
// 3D //
////////

void testMicroRefinement3D() {}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testMicroRefinement2D();
}
