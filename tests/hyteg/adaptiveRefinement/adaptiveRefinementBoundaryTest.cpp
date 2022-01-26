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

using walberla::int64_t;
using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

template < uint_t K, class K_Simplex >
void adaptiveRefinementBoundaryTest( uint_t n_refinements )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Boundary test " << K << "d" );

   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( K == 2 )
   {
      meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   }
   else
   {
      meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
   }

   SetupPrimitiveStorage setupStorage_init( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage_init.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto n_boundary_edges = setupStorage_init.getNumEdgesOnBoundary();
   auto n_boundary_faces = setupStorage_init.getNumFacesOnBoundary();

   adaptiveRefinement::K_Mesh< K_Simplex > mesh( setupStorage_init );

   for ( uint_t ref = 0; ref < n_refinements; ++ref )
   {
      std::vector< PrimitiveID > to_refine;
      for ( auto& el : mesh.get_elements() )
      {
         to_refine.push_back( el->getPrimitiveID() );
      }

      auto& setupStorage = mesh.refineRG( to_refine );

      std::stringstream ss;
      setupStorage.toStream( ss, false );
      WALBERLA_LOG_INFO_ON_ROOT( "Refinement " << ( ref + 1 ) << ": " << ss.str() );

      n_boundary_edges = n_boundary_edges * 2 + n_boundary_faces * 3;
      n_boundary_faces = n_boundary_faces * 4;
      WALBERLA_CHECK_EQUAL( setupStorage.getNumEdgesOnBoundary(), n_boundary_edges );
      WALBERLA_CHECK_EQUAL( setupStorage.getNumFacesOnBoundary(), n_boundary_faces );

      hyteg::SetupPrimitiveStorage::PrimitiveMap primitiveMap;
      setupStorage.getSetupPrimitives( primitiveMap );

      for ( auto& [id, primitive] : primitiveMap )
      {
         auto   flag     = primitive->getMeshBoundaryFlag();
         uint_t boundary = setupStorage.onBoundary( PrimitiveID( id ), true ) ? 1 : 0;
         WALBERLA_CHECK_EQUAL( flag, boundary );
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::adaptiveRefinementBoundaryTest< 2, hyteg::adaptiveRefinement::Simplex2 >( 4 );
   hyteg::adaptiveRefinementBoundaryTest< 3, hyteg::adaptiveRefinement::Simplex3 >( 3 );
}