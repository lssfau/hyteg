/*
 * Copyright (c) 2022 Benjamin Mann
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
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/adaptiverefinement/mesh.hpp"
#include "hyteg/adaptiverefinement/simplex.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"

using walberla::int64_t;
using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

template < uint_t K, class K_Simplex >
void adaptiveRefinementBlendingTest( uint_t n_refinements )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( K == 2 )
   {
      meshInfo = MeshInfo::meshAnnulus( 1, 2, MeshInfo::CRISS, 5, 2 );
   }
   else
   {
      meshInfo = MeshInfo::meshSphericalShell( 3, 2, 1, 2 );
   }

   SetupPrimitiveStorage setupStorage_init( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( K == 2 )
   {
      AnnulusMap::setMap( setupStorage_init );
   }
   else
   {
      IcosahedralShellMap::setMap( setupStorage_init );
   }

   adaptiveRefinement::K_Mesh< K_Simplex > mesh( setupStorage_init );

   for ( uint_t ref = 0; ref < n_refinements; ++ref )
   {
      std::vector< PrimitiveID > to_refine;
      for ( auto& el : mesh.get_elements() )
      {
         to_refine.push_back( el->getPrimitiveID() );
      }

      auto& setupStorage = mesh.refineRG( to_refine );

      auto vertices = setupStorage.getVertices();

      for ( auto& [id, vtx] : vertices )
      {
         Point3D x_blended;
         vtx->getGeometryMap()->evalF( vtx->getCoordinates(), x_blended );
         auto r = x_blended.norm();

         if ( setupStorage.onBoundary( PrimitiveID( id ) ) )
         {
            if ( r < 1.5 )
            {
               WALBERLA_CHECK_FLOAT_EQUAL( r, 1.0 );
            }
            else
            {
               WALBERLA_CHECK_FLOAT_EQUAL( r, 2.0 );
            }
         }
         else
         {
            WALBERLA_CHECK_GREATER( r, 1.0 );
            WALBERLA_CHECK_LESS( r, 2.0 );
         }
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::adaptiveRefinementBlendingTest< 2, hyteg::adaptiveRefinement::Simplex2 >( 2 );
   hyteg::adaptiveRefinementBlendingTest< 3, hyteg::adaptiveRefinement::Simplex3 >( 2 );
}