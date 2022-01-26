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
void adaptiveRefinementTest()
{
   MeshInfo                    meshInfo = MeshInfo::emptyMeshInfo();
   Point3D                     barycenter_0;
   real_t                      volume;
   std::pair< real_t, real_t > min_max_angle = { walberla::math::pi / 4.0, walberla::math::pi / 2.0 };

   uint_t n_el_1, n_el_2, n_el_green_2;

   if ( K == 2 )
   {
      meshInfo     = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
      barycenter_0 = ( 1. / 3. ) * hyteg::Point3D( { 1, 1, 0 } );
      volume       = 1. / 2.;
      n_el_1       = 4;
      n_el_2       = 8;
      n_el_green_2 = 2;
   }
   else
   {
      meshInfo     = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
      barycenter_0 = ( 1. / 4. ) * hyteg::Point3D( { 1, 1, 1 } );
      volume       = 1. / 6.;
      n_el_1       = 8;
      n_el_2       = 20;
      n_el_green_2 = 8;
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   adaptiveRefinement::K_Mesh< K_Simplex > mesh( setupStorage );

   //== test simplex functions before refinement ===

   auto simplex_0  = *( mesh.get_elements().begin() );
   auto vertices_0 = mesh.get_vertices();

   // test Simplex::has_children()
   WALBERLA_CHECK( !simplex_0->has_children() );

   // test Simplex::has_vertex()
   for ( uint_t vtx = 0; vtx < vertices_0.size(); ++vtx )
      WALBERLA_CHECK( simplex_0->has_vertex( vtx ) );
   WALBERLA_CHECK( !simplex_0->has_vertex( 17 ) );

   // test Simplex::barycenter()
   auto barycenter = simplex_0->barycenter( vertices_0 );
   for ( uint_t i = 0; i < 3; ++i )
      WALBERLA_CHECK_FLOAT_EQUAL( barycenter[i], barycenter_0[i] );

   // test Simplex::volume()
   WALBERLA_CHECK_FLOAT_EQUAL( simplex_0->volume( vertices_0 ), volume );

   // test Simplex::min_max_angle()
   auto angles = simplex_0->min_max_angle( vertices_0 );
   WALBERLA_CHECK_FLOAT_EQUAL( min_max_angle.first, angles.first );
   WALBERLA_CHECK_FLOAT_EQUAL( min_max_angle.second, angles.second );

   //=== refine mesh at simplex_0 (the only element) ===

   std::vector< PrimitiveID > to_refine( 1 );
   to_refine[0] = simplex_0->getPrimitiveID();
   mesh.refineRG( to_refine );

   // check number of elements
   WALBERLA_CHECK_EQUAL( mesh.get_elements().size(), n_el_1 );

   // test Simplex::has_children()
   WALBERLA_CHECK( simplex_0->has_children() );

   // test Simplex::has_green_edge()
   for ( auto& el : mesh.get_elements() )
      WALBERLA_CHECK( !el->has_green_edge() );

   //=== refine mesh at simplex_1 (the element sitting at the origin) ===

   uint_t origin = 0; // (0, 0, 0)
   for ( uint_t i = 0; i < vertices_0.size(); ++i )
   {
      if ( vertices_0[i].normSq() <= 0.0 )
      {
         origin = i;
         break;
      }
   }

   std::shared_ptr< K_Simplex > simplex_1 = nullptr; // element at origin
   for ( auto& el : mesh.get_elements() )
   {
      if ( el->has_vertex( origin ) )
      {
         simplex_1 = el;
      }
   }

   to_refine[0] = simplex_1->getPrimitiveID();
   mesh.refineRG( to_refine );

   // check number of elements
   WALBERLA_CHECK_EQUAL( mesh.get_elements().size(), n_el_2 );

   // check number of green elements
   uint_t n_green = 0;
   for ( auto& el : mesh.get_elements() )
   {
      if ( el->has_green_edge() )
         ++n_green;
   }
   WALBERLA_CHECK_EQUAL( n_green, n_el_green_2 );

   // test K_Mesh::volume()
   WALBERLA_CHECK_FLOAT_EQUAL( mesh.volume(), volume );

   // test K_mesh::n_elements()
   WALBERLA_CHECK_EQUAL( mesh.n_elements(), n_el_2 );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::adaptiveRefinementTest< 2, hyteg::adaptiveRefinement::Simplex2 >();
   hyteg::adaptiveRefinementTest< 3, hyteg::adaptiveRefinement::Simplex3 >();
}