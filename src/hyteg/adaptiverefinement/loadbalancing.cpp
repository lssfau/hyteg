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

#include <core/Format.hpp>
#include <core/logging/all.h>
#include <core/mpi/Broadcast.h>
#include <core/mpi/Reduce.h>
#include <numeric>
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

void loadbalancing( const std::vector< Point3D >&      coordinates,
                    std::vector< VertexData >&         vtxs,
                    std::vector< EdgeData >&           edges,
                    std::vector< FaceData >&           faces,
                    std::vector< CellData >&           cells,
                    const std::vector< Neighborhood >& nbrHood,
                    const uint_t&                      n_processes,
                    const uint_t&                      rank )
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

   // unassign everything
   for ( auto& p : vtxs )
   {
      p.setTargetRank( n_processes );
   }
   for ( auto& p : edges )
   {
      p.setTargetRank( n_processes );
   }
   for ( auto& p : faces )
   {
      p.setTargetRank( n_processes );
   }
   for ( auto& p : cells )
   {
      p.setTargetRank( n_processes );
   }

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

      // we only prescribe a maximum for volume elements
      if ( pt == VOL )
      {
         n_max[pt] = n_min + ( ( 0 < mod ) ? 1 : 0 );
      }
      else
      {
         n_max[pt] = n_prim[pt];
      }
   }

   // compute neighboring volume primitives of all primitives
   std::vector< std::vector< uint_t > > nbrVolumes( n_prim[ALL] );
   for ( uint_t idx = 0; idx < nbrHood.size(); ++idx )
   {
      uint_t i = id0[VOL] + idx;
      for ( PT pt : VEFC )
      {
         for ( uint_t j : nbrHood[idx][pt] )
         {
            nbrVolumes[j].push_back( i );
         }
      }
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

   // which primitives are currently assigned to a cluster
   std::vector< bool > isAssigned( n_prim[ALL] + 1, false );
   // how many primitives of each type are assigned to each process
   std::vector< std::array< uint_t, ALL + 1 > > n_assigned( n_processes + 1, std::array< uint_t, ALL + 1 >{} );
   // volume elements assigned to each cluster
   std::vector< std::vector< uint_t > > volume_elements( n_processes );

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
   // compute potential volume of cluster built around element i
   auto predict_volume = [&]( uint_t i ) -> uint_t {
      std::vector< uint_t > Q, Q_new;
      std::vector< bool >   visited( n_prim[ALL], false );
      uint_t                v = 0;
      Q_new.push_back( i );
      visited[i] = true;
      // breadth first search to compute the number of free elements before hitting another cluster
      while ( !Q_new.empty() )
      {
         v += Q_new.size();
         std::swap( Q, Q_new );
         Q_new.clear();

         for ( auto j : Q )
         {
            for ( auto n : nbrVolumes[j] )
            {
               if ( !visited[n] )
               {
                  if ( isAssigned[n] )
                  {
                     return v;
                  }
                  Q_new.push_back( n );
                  visited[n] = true;
               }
            }
         }
      }
      return v;
   };

   // IDs of initial elements
   std::vector< uint_t > initID( n_processes );
   // select initial elements at random
   std::iota( initID.begin(), initID.end(), id0[VOL] );
   std::fill( std::next( isAssigned.begin(), id0[VOL] ), std::next( isAssigned.begin(), id0[VOL] + n_processes ), true );
   // loop over all clusters k and choose initial element maximizing potential cluster size
   for ( uint_t k = 0; k < n_processes; ++k )
   {
      uint_t max_i = 0;           // max_i predict_volume(i)
      uint_t i_max = n_prim[ALL]; // arg max_i predict_volume(i)

      isAssigned[initID[k]] = false;

      // loop over all volume elements
      for ( uint_t i = begin[VOL]; i < end[VOL]; ++i )
      {
         // skip elements that are occupied by another cluster
         if ( isAssigned[i] )
            continue;

         auto v_i = predict_volume( i );

         // find maximum over i
         if ( v_i > max_i )
         {
            max_i = v_i;
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
      initID[k]         = i_max;
      isAssigned[i_max] = true;
   }

   // reset flag
   std::fill( isAssigned.begin(), isAssigned.end(), false );

   // add initial elements
   for ( uint_t k = 0; k < n_processes; ++k )
   {
      assign( initID[k], k );
   }

   // force of attraction between a cluster and its elements
   std::vector< uint_t > attraction( n_processes + 1, n_processes );
   attraction[n_processes] = 0;

   // compute the number of neighbors of i, assigned to k
   auto connectivity = [&]( uint_t i, uint_t k ) -> int {
      if ( k >= n_processes )
      {
         return 0;
      }
      int conn = 0;
      for ( uint_t nbr : nbrVolumes[i] )
      {
         if ( assigned_to( nbr ) == k )
            ++conn;
      }
      return conn;
   };
   // find next primitive of type pt for cluster k
   auto find_next = [&]( uint_t k, PT pt ) -> std::pair< uint_t, int > {
      std::pair< uint_t, int > bestFit = { n_prim[ALL], 0 };

      // loop over all volume elements owned by k
      for ( uint_t id : volume_elements[k] )
      {
         uint_t elIdx = id - id0[VOL];
         // loop over all primitives i of type pt in the neighborhood of el
         for ( uint_t i : nbrHood[elIdx][pt] )
         {
            // i is currently assigned to j
            auto j = assigned_to( i );

            // skip elements already assigned to k
            if ( j == k )
               continue;

            // j binds stronger than k
            if ( attraction[k] < attraction[j] )
               continue;

            // j has only 1 volume element
            if ( j < n_processes && n_assigned[j][VOL] == 1u )
               continue;

            // k is saturated
            if ( n_assigned[k][pt] >= n_max[pt] )
            {
               // i is considered if k binds stronger than j
               if ( attraction[j] < attraction[k] )
                  return { i, 0 };
               else
                  continue;
            }

            // total improvement when assigning i to k
            int improvement = connectivity( i, k ) - connectivity( i, j );
            // bonus due to stronger attraction
            if ( attraction[j] < attraction[k] )
               improvement += 100;

            // i is new best fit
            if ( improvement > bestFit.second )
            {
               bestFit.first  = i;
               bestFit.second = improvement;
            }
         }
      }

      return bestFit;
   };

   // add remaining primitives
   for ( PT pt : CFEV )
   {
      while ( n_assigned[n_processes][pt] < n_prim[pt] )
      {
         // find best fitting element i for cluster k=rank
         auto next_rnk = find_next( rank, pt );

         // collect all elements assigned in the current iteration
         std::map< uint_t, int > bestFit;

         // loop over all clusters k
         for ( uint_t k = 0; k < n_processes; ++k )
         {
            auto next = next_rnk;
            walberla::mpi::broadcastObject( next, int( k ) );
            auto i = next.first;

            // no suitable element found
            if ( primitiveType( i ) != pt )
            {
               continue;
            }

            if ( n_assigned[k][pt] < n_max[pt] ) // cluster k is not saturated
            {
               /* if i has been assigned to another cluster during the same iteration
                     we only reassign it to k if it this is the better fit
                  */
               auto previous = bestFit.find( next.first );
               if ( previous == bestFit.end() || previous->second < next.second )
               {
                  bestFit[next.first] = next.second;
               }
               else
               {
                  continue;
               }

               // unassign element i from its current cluster j
               auto j = unassign( i );
               // j just gave up an element -> restore the forces of attraction
               if ( j < n_processes && pt == VOL )
               {
                  auto threshold = attraction[j];
                  for ( uint_t m = 0; m < n_processes; ++m )
                  {
                     if ( attraction[m] >= threshold )
                        attraction[m] = n_processes;
                  }
               }

               // assign element i to cluster k
               assign( i, k );
            }
            else if ( pt == VOL ) // cluster k is saturated
            {
               /* there are free elements in the nbrhood, but k is saturated already
                  -> k must donate some elements before we can continue.
                  Therefore, we decrease its force of attraction.
               */
               attraction[k] = std::min( attraction[assigned_to( i )] + 1, n_processes );
            }
         }
      }

      /* Those clusters with fewer volume elements should be assigned more
         non-volume elements. We adjust the forces of attraction of each
         cluster accordingly.
      */
      if ( pt == VOL )
      {
         for ( uint_t k = 0; k < n_processes; ++k )
         {
            attraction[k] = 1 + n_max[VOL] - n_assigned[k][VOL];
         }
      }
   }

   WALBERLA_DEBUG_SECTION()
   {
      const uint_t k = rank;

      std::array< uint_t, ALL > n_primitives_k{};

      // count primitives per rank and check that everything got assigned
      for ( auto& p : vtxs )
      {
         WALBERLA_ASSERT_LESS( p.getTargetRank(), n_processes );
         if ( p.getTargetRank() == k )
            ++n_primitives_k[VTX];
      }
      for ( auto& p : edges )
      {
         WALBERLA_ASSERT_LESS( p.getTargetRank(), n_processes );
         if ( p.getTargetRank() == k )
            ++n_primitives_k[EDGE];
      }
      for ( auto& p : faces )
      {
         WALBERLA_ASSERT_LESS( p.getTargetRank(), n_processes );
         if ( p.getTargetRank() == k )
            ++n_primitives_k[FACE];
      }
      for ( auto& p : cells )
      {
         WALBERLA_ASSERT_LESS( p.getTargetRank(), n_processes );
         if ( p.getTargetRank() == k )
            ++n_primitives_k[CELL];
      }

      std::array< uint_t, ALL >      p_min{};
      std::array< uint_t, ALL >      p_max{};
      std::array< real_t, ALL >      p_mean{};
      std::array< std::string, ALL > p_name{ "vertex", "edge", "face", "cell" };

      WALBERLA_LOG_INFO_ON_ROOT( "Primitive distribution:" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s |%10s |%10s |%10s |", "type", "min", "max", "mean" ) );
      for ( PT pt = VTX; pt < ALL; pt = PT( pt + 1 ) )
      {
         p_min[pt]  = walberla::mpi::reduce( n_primitives_k[pt], walberla::mpi::MIN );
         p_max[pt]  = walberla::mpi::reduce( n_primitives_k[pt], walberla::mpi::MAX );
         p_mean[pt] = real_t( walberla::mpi::reduce( n_primitives_k[pt], walberla::mpi::SUM ) ) / real_t( n_processes );
         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( "%20s |%10d |%10d |%10.1f |", p_name[pt].c_str(), p_min[pt], p_max[pt], p_mean[pt] ) );
      }

      // compute surface of each cluster to estimate the qualitiy of the clusters
      uint_t surface_k = 0;
      for ( auto& vol : volume_elements[k] )
      {
         surface_k += ( 1 + VOL );

         for ( auto& nbr : nbrVolumes[vol] )
         {
            if ( assigned_to( nbr ) == k )
               --surface_k;
         }
      }
      auto s_min  = walberla::mpi::reduce( surface_k, walberla::mpi::MIN );
      auto s_max  = walberla::mpi::reduce( surface_k, walberla::mpi::MAX );
      auto s_mean = real_t( walberla::mpi::reduce( surface_k, walberla::mpi::SUM ) ) / real_t( n_processes );

      WALBERLA_LOG_INFO_ON_ROOT( "Quality of clusters:" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s |%10s |%10s |%10s |", "", "min", "max", "mean" ) );
      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%20s |%10d |%10d |%10.1f |", "clustervolume", p_min[VOL], p_max[VOL], p_mean[VOL] ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%20s |%10d |%10d |%10.1f |", "clustersurface", s_min, s_max, s_mean ) );
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
