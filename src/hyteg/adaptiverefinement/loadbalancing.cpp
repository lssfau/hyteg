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

#include <core/logging/all.h>
#include <core/mpi/Broadcast.h>
#include <core/mpi/Reduce.h>
#include <utility>
#include <vector>

#include "simplexData.hpp"

namespace hyteg {
namespace adaptiveRefinement {

/* apply loadbalancing directly on our datastructures */
void loadbalancing( std::vector< VertexData >& vtxs,
                    std::vector< EdgeData >&   edges,
                    std::vector< FaceData >&   faces,
                    std::vector< CellData >&   cells,
                    const uint_t&              n_processes )
{
   // roundrobin
   uint_t i = 0;

   for ( auto& vtx : vtxs )
   {
      vtx.setTargetRank( i % n_processes );
      ++i;
   }
   for ( auto& edge : edges )
   {
      edge.setTargetRank( i % n_processes );
      ++i;
   }
   for ( auto& face : faces )
   {
      face.setTargetRank( i % n_processes );
      ++i;
   }
   for ( auto& cell : cells )
   {
      cell.setTargetRank( i % n_processes );
      ++i;
   }
}

void loadbalancing( const std::vector< Point3D >& coordinates,
                    std::vector< VertexData >&    vtxs,
                    std::vector< EdgeData >&      edges,
                    std::vector< FaceData >&      faces,
                    std::vector< CellData >&      cells,
                    const uint_t&                 n_processes,
                    const uint_t&                 rank )
{
   using PT = PrimitiveType;
   constexpr std::array< PT, ALL > VEFC{ VTX, EDGE, FACE, CELL };
   // constexpr std::array< PT, ALL > CFEV{ CELL, FACE, EDGE, VTX };

   const PT VOL = ( cells.size() == 0 ) ? FACE : CELL;

   // number of primitives of each type
   std::array< uint_t, ALL + 1 > n_prim;
   n_prim[VTX]  = vtxs.size();
   n_prim[EDGE] = edges.size();
   n_prim[FACE] = faces.size();
   n_prim[CELL] = cells.size();
   n_prim[ALL]  = n_prim[VTX] + n_prim[EDGE] + n_prim[FACE] + n_prim[CELL];

   // first Primitive ID for each primitive type
   std::array< uint_t, ALL + 1 > id0{};
   for ( auto pt : VEFC )
   {
      id0[pt + 1] = id0[pt] + n_prim[pt];
   }

   /* We assume that the elements in the input vectors are ordered by
      PrimitiveID and that for each vertex v, edge e, face f and cell c it holds
               id_v < id_e < id_f < id_c
   */
   uint_t check_id = 0;
   auto   check    = [&]( PrimitiveID id ) {
      if ( id.getID() != check_id )
      {
         WALBERLA_ABORT( "Wrong numbering of primitives!" );
      }
      ++check_id;
   };
   for ( auto& p : vtxs )
   {
      check( p.getPrimitiveID() );
   }
   for ( auto& p : edges )
   {
      check( p.getPrimitiveID() );
   }
   for ( auto& p : faces )
   {
      check( p.getPrimitiveID() );
   }
   for ( auto& p : cells )
   {
      check( p.getPrimitiveID() );
   }

   // we only use this algorithm if there are more volume elements than processes
   if ( n_prim[VOL] < n_processes || n_processes < 2 )
   {
      return loadbalancing( vtxs, edges, faces, cells, n_processes );
   }

   // get primitive type of id
   auto primitiveType = [&]( uint_t id ) -> PT {
      PT pt = VTX;
      while ( pt < ALL && id >= id0[pt + 1] )
      {
         pt = PT( pt + 1 );
      }
      return pt;
   };

   // max number of primitives on one rank for each primitive type
   std::array< uint_t, ALL > n_max;
   // distributed id range for each primitive type
   std::array< uint_t, ALL > begin, end;
   for ( auto pt : VEFC )
   {
      auto n_min = n_prim[pt] / n_processes;
      auto mod   = n_prim[pt] % n_processes;

      begin[pt] = id0[pt] + n_min * rank + ( ( rank < mod ) ? rank : mod );
      end[pt]   = begin[pt] + n_min + ( ( rank < mod ) ? 1 : 0 );

      n_max[pt] = n_min + ( ( 0 < mod ) ? 1 : 0 );
   }

   // compute barycenter, radius and volume of all primitives
   std::vector< Point3D > barycenter( n_prim[ALL] );
   std::vector< real_t >  radius( n_prim[ALL] );
   std::vector< real_t >  volume( n_prim[ALL] );
   for ( auto& p : vtxs )
   {
      barycenter[p.getPrimitiveID().getID()] = coordinates[p.getPrimitiveID().getID()];
   }
   for ( auto& p : edges )
   {
      barycenter[p.getPrimitiveID().getID()] = Simplex1::barycenter( p.get_coordinates( coordinates ) );
   }
   for ( auto& p : faces )
   {
      barycenter[p.getPrimitiveID().getID()] = Simplex2::barycenter( p.get_coordinates( coordinates ) );
      radius[p.getPrimitiveID().getID()]     = Simplex2::radius( p.get_coordinates( coordinates ) );
      volume[p.getPrimitiveID().getID()]     = Simplex2::volume( p.get_coordinates( coordinates ) );
   }
   for ( auto& p : cells )
   {
      barycenter[p.getPrimitiveID().getID()] = Simplex3::barycenter( p.get_coordinates( coordinates ) );
      radius[p.getPrimitiveID().getID()]     = Simplex3::radius( p.get_coordinates( coordinates ) );
      volume[p.getPrimitiveID().getID()]     = Simplex3::volume( p.get_coordinates( coordinates ) );
   }

   // barycenter of the primitive cluster corresponding to each process
   std::vector< Point3D > clusterCenter( n_processes );
   // radius of the primitive cluster corresponding to each process
   std::vector< real_t > clusterRadius( n_processes );
   // which primitives are assigned to some process
   std::vector< bool > isAssigned( n_prim[ALL] + 1, false );
   // how many primitives of each type are assigned to each process
   std::vector< std::array< uint_t, ALL + 1 > > n_assigned( n_processes + 1, std::array< uint_t, ALL + 1 >{} );
   // assign primitive i to cluster k
   auto assign = [&]( uint_t i, uint_t k ) -> bool {
      if ( isAssigned[i] )
      {
         return false;
      }

      PT     pt  = primitiveType( i );
      uint_t idx = i - id0[pt];
      if ( pt == VTX )
      {
         vtxs[idx].setTargetRank( k );
      }
      else if ( pt == EDGE )
      {
         edges[idx].setTargetRank( k );
      }
      else if ( pt == FACE )
      {
         faces[idx].setTargetRank( k );
      }
      else if ( pt == CELL )
      {
         cells[idx].setTargetRank( k );
      }
      else
      {
         return false;
      }
      // mark as assigned
      ++n_assigned[k][pt];
      ++n_assigned[k][ALL];
      ++n_assigned[n_processes][pt];
      ++n_assigned[n_processes][ALL];
      isAssigned[i] = true;

      // update center and radius
      if ( pt == VOL )
      {
         clusterCenter[k] *= real_t( n_assigned[k][VOL] - 1 ); // undo previous scaling
         clusterCenter[k] += barycenter[i];                    // add barycenter of new element
         clusterCenter[k] /= real_t( n_assigned[k][VOL] );     // apply scaling

         clusterRadius[k] = std::max( clusterRadius[k], ( clusterCenter[k] - barycenter[i] ).norm() + radius[i] );
      }

      return true;
   };

   // unassign primitive i from its current cluster
   auto unassign = [&]( uint_t i ) -> uint_t {
      if ( !isAssigned[i] )
      {
         return n_processes;
      }

      PT     pt  = primitiveType( i );
      uint_t idx = i - id0[pt];
      uint_t k   = n_processes;
      if ( pt == VTX )
      {
         k = vtxs[idx].getTargetRank();
         vtxs[idx].setTargetRank( n_processes );
      }
      else if ( pt == EDGE )
      {
         k = edges[idx].getTargetRank();
         edges[idx].setTargetRank( n_processes );
      }
      else if ( pt == FACE )
      {
         k = faces[idx].getTargetRank();
         faces[idx].setTargetRank( n_processes );
      }
      else if ( pt == CELL )
      {
         k = cells[idx].getTargetRank();
         cells[idx].setTargetRank( n_processes );
      }
      else
      {
         return n_processes;
      }

      // mark as unassigned
      --n_assigned[k][pt];
      --n_assigned[k][ALL];
      --n_assigned[n_processes][pt];
      --n_assigned[n_processes][ALL];
      isAssigned[i] = false;

      // update center
      if ( pt == VOL )
      {
         clusterCenter[k] *= real_t( n_assigned[k][VOL] + 1 ); // undo previous scaling
         clusterCenter[k] -= barycenter[i];                    // substract barycenter of removed element
         clusterCenter[k] /= real_t( n_assigned[k][VOL] );     // apply scaling
      }

      return k;
   };

   // which rank is primitive i currently assigned to
   auto assigned_to = [&]( uint_t i ) -> uint_t {
      if ( !isAssigned[i] )
      {
         return n_processes;
      }

      PT     pt  = primitiveType( i );
      uint_t idx = i - id0[pt];
      if ( pt == VTX )
      {
         return vtxs[idx].getTargetRank();
      }
      else if ( pt == EDGE )
      {
         return edges[idx].getTargetRank();
      }
      else if ( pt == FACE )
      {
         return faces[idx].getTargetRank();
      }
      else if ( pt == CELL )
      {
         return cells[idx].getTargetRank();
      }
      else
      {
         return n_processes;
      }
   };

   // IDs of initial elements
   std::vector< uint_t > initID( n_processes, n_prim[ALL] );

   /* select initial elements for each cluster to maximize
       min_j!=k||clusterCenter[j] - clusterCenter[k]||
   */
   for ( uint_t cycle = 0; cycle < 100; ++cycle ) // backup stopping criterion in case of cyclic states
   {
      bool stateChange = false;

      /* loop over all clusters k and move their ceinter to barycenter[i_max]
         for an element i_max maximizing
            min_j ||clusterCenter[j] - barycenter[i]||
      */
      for ( uint_t k = 0; k < n_processes; ++k )
      {
         real_t max_i = 0.0;         // max_i min_j ||clusterCenter[j] - barycenter[i]||
         uint_t i_max = n_prim[ALL]; // arg max_i min_j ||clusterCenter[j] - barycenter[i]||

         // loop over all volume elements
         for ( uint_t i = begin[VOL]; i < end[VOL]; ++i )
         {
            // skip elements that are occupied by another cluster
            if ( isAssigned[i] && i != initID[k] )
               continue;

            // min_j ||clusterCenter[j] - barycenter[i]||
            real_t min_j = std::numeric_limits< real_t >::max();

            // loop over all clusters j
            for ( uint_t j = 0; j < n_processes; ++j )
            {
               // skip j==k and j uninitialized
               if ( j == k || initID[j] == n_prim[ALL] )
               {
                  continue;
               }

               // ||clusterCenter[j] - barycenter[i]||
               auto d = ( barycenter[initID[j]] - barycenter[i] ).norm();

               // find minimum over j
               if ( d < min_j )
               {
                  min_j = d;
               }
            }

            // find maximum over i
            if ( min_j > max_i )
            {
               max_i = min_j;
               i_max = i;
            }
         }
         // global maximum
         auto global_max = walberla::mpi::allReduce( max_i, walberla::mpi::MAX );
         if ( global_max > max_i )
         {
            i_max = n_prim[ALL];
         }
         i_max = walberla::mpi::allReduce( i_max, walberla::mpi::MIN );

         WALBERLA_ASSERT( i_max < n_prim[ALL] );

         // apply changes
         if ( i_max != initID[k] )
         {
            stateChange = true;
            // update cluster k
            isAssigned[initID[k]] = false;
            initID[k]             = i_max;
            isAssigned[initID[k]] = true;
         }
      }

      if ( !stateChange )
      {
         break;
      }
   }

   // reset flag
   std::fill( isAssigned.begin(), isAssigned.end(), false );

   // add initial elements
   for ( uint_t k = 0; k < n_processes; ++k )
   {
      assign( initID[k], k );
   }

   // add remaining volume elements
   std::vector< std::vector< uint_t > > stole_from( n_processes );
   while ( n_assigned[n_processes][VOL] < n_prim[VOL] )
   {
      bool activity = false;
      // loop over all clusters k
      for ( uint_t k = 0; k < n_processes; ++k )
      {
         if ( n_assigned[k][VOL] >= n_max[VOL] )
            continue;

         // min_i ||clusterCenter[k] - barycenter[i]||
         real_t min_i[2] = { clusterRadius[k], clusterRadius[k] };
         // arg min_i ||clusterCenter[k] - barycenter[i]||
         uint_t i_min[2] = { n_prim[ALL], n_prim[ALL] };

         // loop over all volume elements i to find minimum
         for ( uint_t i = begin[VOL]; i < end[VOL]; ++i )
         {
            auto j = assigned_to( i );

            // skip elements already assigned to k
            if ( j == k )
               continue;

            if ( std::find( stole_from[k].begin(), stole_from[k].end(), j ) != stole_from[k].end() )
               continue;

            // ||clusterCenter[k] - barycenter[i]||
            auto d_ik = ( clusterCenter[k] - barycenter[i] ).norm();

            // prioritize
            auto prio = uint_t( isAssigned[i] );

            // find minimum over i
            if ( d_ik < min_i[prio] &&                                   // d_ik is new minimum AND
                 ( prio == 0u                                            // (   i is still free OR
                   || d_ik < ( clusterCenter[j] - barycenter[i] ).norm() //     i is closer to k than to j OR
                   || n_assigned[j][VOL] > n_assigned[k][VOL] ) )        //     j has more elements than k )
            //  || n_assigned[j][VOL] == n_max[VOL] ) )        //     j has more elements than k )
            {
               min_i[prio] = d_ik;
               i_min[prio] = i;
            }
         }

         // global minimum
         real_t min_i_global[2];
         for ( auto prio : { 0u, 1u } )
         {
            min_i_global[prio] = walberla::mpi::allReduce( min_i[prio], walberla::mpi::MIN );
         }

         auto which   = ( min_i_global[0] < clusterRadius[k] ) ? 0u : 1u;
         i_min[which] = ( min_i_global[which] < min_i[which] ) ? n_prim[ALL] : i_min[which];
         i_min[which] = walberla::mpi::allReduce( i_min[which], walberla::mpi::MIN );

         if ( i_min[which] < id0[VOL + 1] )
         {
            // unassign element from its current cluster if necessary
            auto j = unassign( i_min[which] );
            if ( j < n_processes )
            {
               stole_from[j].push_back( k );
            }
            // assign element i_min to cluster k_min
            assign( i_min[which], k );

            WALBERLA_LOG_INFO_ON_ROOT( "steal: " << j << " --" << i_min[which] << "--> " << k );
            activity = true;
         }

         // stole_from[k].clear();
      }

      auto min = walberla::mpi::allReduce( n_assigned[rank][VOL], walberla::mpi::MIN );
      auto max = walberla::mpi::allReduce( n_assigned[rank][VOL], walberla::mpi::MAX );

      WALBERLA_LOG_INFO_ON_ROOT( "assigned: " << n_assigned[n_processes][VOL] << "/" << n_prim[VOL] );
      WALBERLA_LOG_INFO_ON_ROOT( "min, max, should: " << min << ", " << max << ", " << n_max[VOL] );

      if ( !activity )
      {
         break;
      }
   }

   // add interface primitives
   for ( auto pt : VEFC )
   {
      if ( n_assigned[n_processes][pt] >= n_prim[pt] )
      {
         continue;
      }

      // todo swap loops

      // "loop" over all clusters k (MPI)
      uint_t k = rank;

      // loop over all elements i of type pt
      for ( uint_t i = id0[pt]; i < id0[pt + 1]; ++i )
      {
         // skip assigned elements
         if ( isAssigned[i] )
            continue;

         // ||clusterCenter[k] - barycenter[i]||
         real_t d_ik = std::numeric_limits< real_t >::max();
         if ( n_assigned[k][pt] < n_max[pt] )
         {
            d_ik = ( clusterCenter[k] - barycenter[i] ).norm();
         }

         // minimum over k
         auto min_k = walberla::mpi::allReduce( d_ik, walberla::mpi::MIN );
         auto k_min = ( min_k < d_ik ) ? n_processes : k;
         k_min      = walberla::mpi::allReduce( k_min, walberla::mpi::MIN );

         // assign element i to cluster k_min
         assign( i, k_min );
      }
   }
}

void loadbalancing( const std::vector< Point3D >& coordinates,
                    std::vector< VertexData >&    vtxs,
                    std::vector< EdgeData >&      edges,
                    std::vector< FaceData >&      faces,
                    std::vector< CellData >&      cells,
                    std::vector< Neighborhood >   nbrHood,
                    const uint_t&                 n_processes,
                    const uint_t&                 rank )
{
   using PT = PrimitiveType;
   constexpr std::array< PT, ALL > VEFC{ VTX, EDGE, FACE, CELL };
   constexpr std::array< PT, ALL > CFEV{ CELL, FACE, EDGE, VTX };

   const PT VOL = ( cells.size() == 0 ) ? FACE : CELL;

   // number of primitives of each type
   std::array< uint_t, ALL + 1 > n_prim;
   n_prim[VTX]  = vtxs.size();
   n_prim[EDGE] = edges.size();
   n_prim[FACE] = faces.size();
   n_prim[CELL] = cells.size();
   n_prim[ALL]  = n_prim[VTX] + n_prim[EDGE] + n_prim[FACE] + n_prim[CELL];

   // first Primitive ID for each primitive type
   std::array< uint_t, ALL + 1 > id0{};
   for ( auto pt : VEFC )
   {
      id0[pt + 1] = id0[pt] + n_prim[pt];
   }

   /* We assume that the elements in the input vectors are ordered by
      PrimitiveID and that for each vertex v, edge e, face f and cell c it holds
               id_v < id_e < id_f < id_c
   */
   uint_t check_id = 0;
   auto   check    = [&]( PrimitiveID id ) {
      if ( id.getID() != check_id )
      {
         WALBERLA_ABORT( "Wrong numbering of primitives!" );
      }
      ++check_id;
   };
   for ( auto& p : vtxs )
   {
      check( p.getPrimitiveID() );
   }
   for ( auto& p : edges )
   {
      check( p.getPrimitiveID() );
   }
   for ( auto& p : faces )
   {
      check( p.getPrimitiveID() );
   }
   for ( auto& p : cells )
   {
      check( p.getPrimitiveID() );
   }

   // we only use this algorithm if there are more volume elements than processes
   if ( n_prim[VOL] < n_processes || n_processes < 2 )
   {
      return loadbalancing( vtxs, edges, faces, cells, n_processes );
   }

   // get primitive type of id
   auto primitiveType = [&]( uint_t id ) -> PT {
      PT pt = VTX;
      while ( pt < ALL && id >= id0[pt + 1] )
      {
         pt = PT( pt + 1 );
      }
      return pt;
   };

   // max number of primitives on one rank for each primitive type
   std::array< uint_t, ALL > n_max;
   // todo: check if this is still required
   // distributed id range for each primitive type
   std::array< uint_t, ALL > begin, end;
   for ( auto pt : VEFC )
   {
      auto n_min = n_prim[pt] / n_processes;
      auto mod   = n_prim[pt] % n_processes;

      begin[pt] = id0[pt] + n_min * rank + ( ( rank < mod ) ? rank : mod );
      end[pt]   = begin[pt] + n_min + ( ( rank < mod ) ? 1 : 0 );

      n_max[pt] = n_min + ( ( 0 < mod ) ? 1 : 0 );
   }

   // todo: check if this is still required
   // compute barycenter, radius and volume of all primitives
   std::vector< Point3D > barycenter( n_prim[ALL] );
   std::vector< real_t >  radius( n_prim[ALL] );
   std::vector< real_t >  volume( n_prim[ALL] );
   for ( auto& p : vtxs )
   {
      barycenter[p.getPrimitiveID().getID()] = coordinates[p.getPrimitiveID().getID()];
   }
   for ( auto& p : edges )
   {
      barycenter[p.getPrimitiveID().getID()] = Simplex1::barycenter( p.get_coordinates( coordinates ) );
   }
   for ( auto& p : faces )
   {
      barycenter[p.getPrimitiveID().getID()] = Simplex2::barycenter( p.get_coordinates( coordinates ) );
      radius[p.getPrimitiveID().getID()]     = Simplex2::radius( p.get_coordinates( coordinates ) );
      volume[p.getPrimitiveID().getID()]     = Simplex2::volume( p.get_coordinates( coordinates ) );
   }
   for ( auto& p : cells )
   {
      barycenter[p.getPrimitiveID().getID()] = Simplex3::barycenter( p.get_coordinates( coordinates ) );
      radius[p.getPrimitiveID().getID()]     = Simplex3::radius( p.get_coordinates( coordinates ) );
      volume[p.getPrimitiveID().getID()]     = Simplex3::volume( p.get_coordinates( coordinates ) );
   }

   std::vector< bool > isAssigned( n_prim[ALL] + 1, false );
   // how many primitives of each type are assigned to each process
   std::vector< std::array< uint_t, ALL + 1 > > n_assigned( n_processes + 1, std::array< uint_t, ALL + 1 >{} );
   // volume elements assigned to each cluster
   std::vector< std::vector< uint_t > > volume_elements( n_processes );
   std::vector< std::vector< uint_t > > history( n_processes );
   // assign primitive i to cluster k
   auto assign = [&]( uint_t i, uint_t k ) -> bool {
      if ( isAssigned[i] )
      {
         return false;
      }

      PT     pt  = primitiveType( i );
      uint_t idx = i - id0[pt];
      if ( pt == VTX )
      {
         vtxs[idx].setTargetRank( k );
      }
      else if ( pt == EDGE )
      {
         edges[idx].setTargetRank( k );
      }
      else if ( pt == FACE )
      {
         faces[idx].setTargetRank( k );
      }
      else if ( pt == CELL )
      {
         cells[idx].setTargetRank( k );
      }
      else
      {
         return false;
      }
      // mark as assigned
      ++n_assigned[k][pt];
      ++n_assigned[k][ALL];
      ++n_assigned[n_processes][pt];
      ++n_assigned[n_processes][ALL];
      isAssigned[i] = true;

      if ( pt == VOL )
      {
         volume_elements[k].push_back( i );
      }

      return true;
   };
   // unassign primitive i from its current cluster
   auto unassign = [&]( uint_t i ) -> uint_t {
      if ( !isAssigned[i] )
      {
         return n_processes;
      }

      PT     pt  = primitiveType( i );
      uint_t idx = i - id0[pt];
      uint_t k   = n_processes;
      if ( pt == VTX )
      {
         k = vtxs[idx].getTargetRank();
         vtxs[idx].setTargetRank( n_processes );
      }
      else if ( pt == EDGE )
      {
         k = edges[idx].getTargetRank();
         edges[idx].setTargetRank( n_processes );
      }
      else if ( pt == FACE )
      {
         k = faces[idx].getTargetRank();
         faces[idx].setTargetRank( n_processes );
      }
      else if ( pt == CELL )
      {
         k = cells[idx].getTargetRank();
         cells[idx].setTargetRank( n_processes );
      }
      if ( k == n_processes )
      {
         return n_processes;
      }

      // mark as unassigned
      --n_assigned[k][pt];
      --n_assigned[k][ALL];
      --n_assigned[n_processes][pt];
      --n_assigned[n_processes][ALL];
      isAssigned[i] = false;

      if ( pt == VOL )
      {
         volume_elements[k].erase( std::find( volume_elements[k].begin(), volume_elements[k].end(), i ) );
         history[k].push_back( i );
      }

      return k;
   };
   // which rank is primitive i currently assigned to
   auto assigned_to = [&]( uint_t i ) -> uint_t {
      if ( !isAssigned[i] )
      {
         return n_processes;
      }

      PT     pt  = primitiveType( i );
      uint_t idx = i - id0[pt];
      if ( pt == VTX )
      {
         return vtxs[idx].getTargetRank();
      }
      else if ( pt == EDGE )
      {
         return edges[idx].getTargetRank();
      }
      else if ( pt == FACE )
      {
         return faces[idx].getTargetRank();
      }
      else if ( pt == CELL )
      {
         return cells[idx].getTargetRank();
      }
      else
      {
         return n_processes;
      }
   };
   // IDs of initial elements
   std::vector< uint_t > initID( n_processes, n_prim[ALL] );
   WALBERLA_LOG_INFO_ON_ROOT("initial elements");
   /* select initial elements for each cluster to maximize
       min_j!=k||clusterCenter[j] - clusterCenter[k]||
   */ // todo: use weighted norm s.th. smaller clusters elements increase the effective distance
   for ( uint_t cycle = 0; cycle < 100; ++cycle ) // backup stopping criterion in case of cyclic states
   {
      bool stateChange = false;

      /* loop over all clusters k and move their ceinter to barycenter[i_max] with
            i_max = arg max_i min_j ||clusterCenter[j] - barycenter[i]||
      */
      for ( uint_t k = 0; k < n_processes; ++k )
      {
         real_t max_i = 0.0;         // max_i min_j ||clusterCenter[j] - barycenter[i]||
         uint_t i_max = n_prim[ALL]; // arg max_i min_j ||clusterCenter[j] - barycenter[i]||

         // loop over all volume elements
         for ( uint_t i = begin[VOL]; i < end[VOL]; ++i )
         {
            // skip elements that are occupied by another cluster
            if ( isAssigned[i] && i != initID[k] )
               continue;

            // min_j ||clusterCenter[j] - barycenter[i]||
            real_t min_j = std::numeric_limits< real_t >::max();

            // loop over all clusters j
            for ( uint_t j = 0; j < n_processes; ++j )
            {
               // skip j==k and j uninitialized
               if ( j == k || initID[j] == n_prim[ALL] )
               {
                  continue;
               }

               // ||clusterCenter[j] - barycenter[i]||
               auto d = ( barycenter[initID[j]] - barycenter[i] ).norm();

               // find minimum over j
               if ( d < min_j )
               {
                  min_j = d;
               }
            }

            // find maximum over i
            if ( min_j > max_i )
            {
               max_i = min_j;
               i_max = i;
            }
         }
         // global maximum
         auto global_max = walberla::mpi::allReduce( max_i, walberla::mpi::MAX );
         if ( global_max > max_i )
         {
            i_max = n_prim[ALL];
         }
         i_max = walberla::mpi::allReduce( i_max, walberla::mpi::MIN );

         WALBERLA_ASSERT( i_max < n_prim[ALL] );

         // apply changes
         if ( i_max != initID[k] )
         {
            stateChange = true;
            // update cluster k
            isAssigned[initID[k]] = false;
            initID[k]             = i_max;
            isAssigned[initID[k]] = true;
         }
      }

      if ( !stateChange )
      {
         break;
      }
   }

   // reset flag
   std::fill( isAssigned.begin(), isAssigned.end(), false );

   // add initial elements
   WALBERLA_LOG_INFO_ON_ROOT("assign");
   for ( uint_t k = 0; k < n_processes; ++k )
   {
      assign( initID[k], k );
   }

   // recall whether a cluster stole from another
   std::vector< std::vector< uint_t > > stole_from( n_processes );
   // find next primitive of type pt for cluster k
   auto find_next = [&]( uint_t k, PT pt ) -> uint_t {
      uint_t steal = n_prim[ALL];

      // loop over all volume elements owned by k
      for ( uint_t id : volume_elements[k] )
      {
         uint_t elIdx = id - id0[VOL];
         // loop over all neighbor primitives i of el
         for ( uint_t i : nbrHood[elIdx][pt] )
         {
            // if i is not owned by any cluster, assign it to k
            if ( !isAssigned[i] )
            {
               return i;
            }

            // i is currently assigned to j
            auto j = assigned_to( i );

            // skip elements already assigned to k
            if ( j == k )
               continue;

            // only steal from bigger clusters
            if ( n_assigned[k][pt] >= n_assigned[j][pt] )
               continue;

            // we do not allow stealing back
            if ( pt == VOL )
            {
               if ( std::find( history[k].begin(), history[k].end(), i ) != history[k].end() )
                  continue;
               if ( std::find( stole_from[k].begin(), stole_from[k].end(), j ) != stole_from[k].end() )
                  continue;
            }

            // mark the first possible element for stealing
            if ( steal == n_prim[ALL] )
               steal = i;
         }
      }

      // if no free element was found, steal the best fitting one
      return steal;
   };

   // add remaining primitives
   for ( PT pt : CFEV )
   {
      WALBERLA_LOG_INFO_ON_ROOT("add elements of type " << pt);

      while ( n_assigned[n_processes][pt] < n_prim[pt] )
      {
         // loop over all clusters k
         for ( uint_t k = 0; k < n_processes; ++k )
         {
            if ( n_assigned[k][pt] >= n_max[pt] )
               continue;

            auto i = find_next( k, pt );

            if ( i < id0[pt + 1] )
            {
               // unassign element from its current cluster if necessary
               auto j = unassign( i );
               if ( j < n_processes )
               {
                  stole_from[j].push_back( k );
               }
               // assign element i to cluster k
               assign( i, k );
            }
         }
         WALBERLA_LOG_INFO_ON_ROOT("assigned: " <<  n_assigned[n_processes][pt] << "/" << n_prim[pt]);

      }
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
