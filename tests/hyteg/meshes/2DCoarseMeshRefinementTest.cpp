/*
* Copyright (c) 2024 Andreas Burkhart
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

#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   uint_t nMacroElements = 1;
   real_t hMax           = 1.0;
   bool   dp             = std::is_same< real_t, double >();

   for ( uint_t i = 0; i < 3; i++ )
   {
      // create mesh
      MeshInfo meshInfo = MeshInfo::refinedCoarseMesh(
          MeshInfo::meshAnnulus( real_c( 1.0 ), real_c( 2.0 ), MeshInfo::meshFlavour::CRISS, 5, 3 ), i );

      SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

      // set load balancing
      loadbalancing::roundRobinVolume( setupStorage );

      // create storage
      auto storage_ = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

      uint_t nMacroElementsPrevious = nMacroElements;
      real_t hMaxPrevious           = hMax;

      if ( storage_->hasGlobalCells() )
      {
         nMacroElements = storage_->getNumberOfGlobalCells();
      }
      else
      {
         nMacroElements = storage_->getNumberOfGlobalFaces();
      }
      hMax = MeshQuality::getMaximalEdgeLength( storage_, 0 );

      WALBERLA_LOG_INFO_ON_ROOT( "---------------------" << i << "---------------------" );
      WALBERLA_LOG_INFO_ON_ROOT( "Macro Elements: " << nMacroElements );
      WALBERLA_LOG_INFO_ON_ROOT( "hMax: " << hMax );

      if ( i > 0 )
      {
         real_t ratio_hMax           = hMax / hMaxPrevious;
         real_t ratio_nMacroElements = real_c( nMacroElements ) / real_c( nMacroElementsPrevious );

         WALBERLA_LOG_INFO_ON_ROOT( "ratio_hMax: " << ratio_hMax );
         WALBERLA_LOG_INFO_ON_ROOT( "ratio_nMacroElements: " << ratio_nMacroElements );

         WALBERLA_CHECK_LESS( std::abs( ratio_hMax - real_c( 0.5 ) ), dp ? 1e-15 : 1e-7 )
         WALBERLA_CHECK_LESS( std::abs( ratio_nMacroElements - real_c( 4.0 ) ), dp ? 1e-15 : 1e-7 )
      }
   }

   return EXIT_SUCCESS;
}