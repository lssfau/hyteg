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
#include <core/mpi/Reduce.h>

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
   // type of primitive
   enum PT
   {
      VTX  = 0,
      EDGE = 1,
      FACE = 2,
      CELL = 3,
      ALL  = 4
   };
   constexpr std::array< PT, ALL > VEFC{ VTX, EDGE, FACE, CELL };
   constexpr std::array< PT, ALL > CFEV{ CELL, FACE, EDGE, VTX };

   const PT VOL = ( cells.size() == 0 ) ? FACE : CELL;

   /* we assume that the elements in the input vectors are ordered by PrimitiveID
      and that for each vertex v, edge e, face f and cell c it holds
               id_v < id_e < id_f < id_c
   */
   std::array< uint_t, ALL + 2 > id0;
   id0[VTX]     = 0;                        // vtxID0
   id0[EDGE]    = id0[VTX] + vtxs.size();   // edgeID0
   id0[FACE]    = id0[EDGE] + edges.size(); // faceID0
   id0[CELL]    = id0[FACE] + faces.size(); // cellID0
   id0[ALL]     = id0[CELL] + cells.size(); // n_all
   id0[ALL + 1] = id0[ALL] + 1;             // n_all+1

   // we only use this algorithm if there are more volume elements than processes
   if ( id0[VOL + 1] - id0[VOL] < n_processes || n_processes < 2 )
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
      auto n_tot = id0[pt + 1] - id0[pt];
      auto n_min = n_tot / n_processes;
      auto mod   = n_tot % n_processes;

      begin[pt] = id0[pt] + n_min * rank + ( ( rank < mod ) ? rank : mod );
      end[pt]   = begin[pt] + n_min + ( ( rank < mod ) ? 1 : 0 );

      n_max[pt] = n_min + ( ( 0 < mod ) ? 1 : 0 );
   }

   // compute barycenter of all primitives
   std::vector< Point3D > barycenter( id0[ALL] + 1 );
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
   }
   for ( auto& p : cells )
   {
      barycenter[p.getPrimitiveID().getID()] = Simplex3::barycenter( p.get_coordinates( coordinates ) );
   }

   // barycenter of the primitive cluster corresponding to each process
   std::vector< Point3D > clusterCenter( n_processes );
   // set value of uninitialized center to max
   barycenter[id0[ALL]].setAll( std::numeric_limits< real_t >::max() );
   // which primitives are assigned to some process
   std::vector< bool > isAssigned( id0[ALL] + 1, false );
   // how many primitives of each type are assigned to each process
   std::vector< std::array< uint_t, ALL + 1 > > n_assigned( n_processes, std::array< uint_t, ALL + 1 >{} );
   // total number of assigned primitives
   uint_t n_assigned_total = 0;
   // assign primitive i to cluster k
   auto assign = [&]( uint_t i, uint_t k ) -> bool {
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
      ++n_assigned_total;
      isAssigned[i] = true;

      // update center
      if ( pt == VOL )
      {
         clusterCenter[k] *= real_t( n_assigned[k][ALL] - 1 ); // undo previous scaling
         clusterCenter[k] += barycenter[i];                    // add barycenter of new point
         clusterCenter[k] /= real_t( n_assigned[k][ALL] );     // apply scaling
      }

      return true;
   };

   // IDs of initial elements
   std::vector< uint_t > initID( n_processes, id0[ALL] );

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
         real_t max_i = 0.0;      // max_i min_j ||clusterCenter[j] - barycenter[i]||
         uint_t i_max = id0[ALL]; // arg max_i min_j ||clusterCenter[j] - barycenter[i]||

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
               if ( j == k || initID[j] == id0[ALL] )
               {
                  continue;
               }

               // ||clusterCenter[j] - barycenter[i]||
               auto d = ( barycenter[initID[j]] - barycenter[i] ).normSq();

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
            i_max = id0[ALL];
         }
         i_max = walberla::mpi::allReduce( i_max, walberla::mpi::MIN );

         WALBERLA_ASSERT( i_max < id0[ALL] );

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

   // add initial elements
   for ( uint_t clusterID = 0; clusterID < n_processes; ++clusterID )
   {
      assign( initID[clusterID], clusterID );
   }

   // add elements
   while ( n_assigned_total < id0[ALL] )
   {
      for ( uint_t clusterID = 0; clusterID < n_processes; ++clusterID )
      {
         // select next primitive
         auto next = id0[ALL];
         // choose element with its barycenter closest to center of cluster
         real_t d_min = std::numeric_limits< real_t >::max();

         for ( auto pt : VEFC )
         {
            if ( n_assigned[clusterID][pt] >= n_max[pt] )
               continue;

            for ( uint_t i = begin[pt]; i < end[pt]; ++i )
            {
               if ( isAssigned[i] )
                  continue;

               auto d = ( clusterCenter[clusterID] - barycenter[i] ).normSq();
               if ( d < d_min )
               {
                  d_min = d;
                  next  = i;
               }
            }
         }
         auto global_min = walberla::mpi::allReduce( d_min, walberla::mpi::MIN );
         if ( global_min < d_min )
         {
            next = id0[ALL];
         }

         next = walberla::mpi::allReduce( next, walberla::mpi::MIN );

         assign( next, clusterID );
      }
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
