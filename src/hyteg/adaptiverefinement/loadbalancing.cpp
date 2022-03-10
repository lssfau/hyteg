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
   for ( auto pt = VTX; pt <= CELL; pt = PT( pt + 1 ) )
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
   std::vector< bool > isAssigned( id0[ALL], false );
   // how many primitives of each type are assigned to each process
   std::vector< std::array< uint_t, ALL + 1 > > n_assigned( n_processes, std::array< uint_t, ALL + 1 >{} );
   // total number of assigned primitives
   uint_t n_assigned_total = 0;
   // assign primitive to cluster
   auto assign = [&]( uint_t id, uint_t clusterID ) -> bool {
      PT pt = primitiveType( id );
      if ( pt == VTX )
      {
         vtxs[id].setTargetRank( clusterID );
      }
      else if ( pt == EDGE )
      {
         edges[id - id0[EDGE]].setTargetRank( clusterID );
      }
      else if ( pt == FACE )
      {
         faces[id - id0[FACE]].setTargetRank( clusterID );
      }
      else if ( pt == CELL )
      {
         cells[id - id0[CELL]].setTargetRank( clusterID );
      }
      else
      {
         return false;
      }
      // mark as assigned
      ++n_assigned[clusterID][pt];
      ++n_assigned[clusterID][ALL];
      ++n_assigned_total;
      isAssigned[id] = true;

      return true;
   };

   // primitive types considered as initial elements
   const PT pt_INIT = primitiveType( n_processes );
   // IDs of initial elements
   std::vector< uint_t > initID( n_processes, id0[ALL] );

   // select initial elements s.th. the distance between the clusters is maximized
   real_t d_sum     = 0.0;
   real_t d_sum_old = -1.0;
   while ( d_sum > d_sum_old )
   {
      d_sum_old = d_sum;

      WALBERLA_LOG_INFO_ON_ROOT( "update initial cluster elements" );
      // update initial element of each cluster
      for ( uint_t clusterID = 0; clusterID < n_processes; ++clusterID )
      {
         /* the initial element of each cluster is chosen s.th.
            id = arg max_i min_j ||barycenter[i] - clusterCenter[j]||
         */
         auto   id    = initID[clusterID];
         real_t d_max = -0.0;

         for ( auto pt = VTX; pt <= pt_INIT; pt = PT( pt + 1 ) )
         {
            for ( uint_t i = begin[pt]; i < end[pt]; ++i )
            {
               // skip elements that are occupied by another cluster
               if ( isAssigned[i] && i != initID[clusterID] )
                  continue;

               // find d_i = min_j ||clusterCenter[j] - barycenter[i]||
               real_t d_i = std::numeric_limits< real_t >::max();
               for ( uint_t j = 0; j < n_processes; ++j )
               {
                  if ( j == clusterID )
                     continue;

                  auto d = ( barycenter[initID[j]] - barycenter[i] ).normSq();
                  if ( d < d_i )
                  {
                     d_i = d;
                  }
               }
               if ( d_i > d_max )
               {
                  d_max = d_i;
                  id    = i;
               }
            }
         }
         auto global_max = walberla::mpi::allReduce( d_max, walberla::mpi::MAX );
         if ( global_max > d_max )
         {
            id = id0[ALL];
         }

         id = walberla::mpi::allReduce( id, walberla::mpi::MIN );

         if ( id < id0[ALL] && id != initID[clusterID] )
         {
            auto initID_old   = initID[clusterID];
            initID[clusterID] = id;

            // check if this is actually an improvement
            real_t d_sum_tst = 0.0;
            for ( uint_t k = 0; k < n_processes; ++k )
            {
               // find d_k = min_j ||clusterCenter[j] - clusterCenter[k]|| for j!=k
               real_t d_k = std::numeric_limits< real_t >::max();
               for ( uint_t j = 0; j < n_processes; ++j )
               {
                  if ( j == k )
                     continue;

                  auto d = ( barycenter[initID[j]] - barycenter[initID[k]] ).normSq();
                  if ( d < d_k )
                  {
                     d_k = d;
                  }
               }

               d_sum_tst += d_k;
            }

            if ( d_sum_tst > d_sum_old )
            {
               d_sum = d_sum_tst;
               WALBERLA_LOG_INFO_ON_ROOT( "clusterID: " << clusterID << " updated from " << initID_old << " to " << id );

               // swap initial element for this one
               if ( initID_old < id0[ALL] )
               {
                  isAssigned[initID_old] = false;
               }
               isAssigned[id] = true;
            }
            else
            {
               initID[clusterID] = initID_old;
            }
         }
      }

      WALBERLA_LOG_INFO_ON_ROOT( "sum(d) = " << d_sum );
   }

   // add initial elements
   for ( uint_t clusterID = 0; clusterID < n_processes; ++clusterID )
   {
      if ( assign( initID[clusterID], clusterID ) )
      {
         clusterCenter[clusterID] = barycenter[initID[clusterID]];
      }
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

         for ( auto pt = VTX; pt <= CELL; pt = PT( pt + 1 ) )
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

         if ( assign( next, clusterID ) )
         {
            // update center
            clusterCenter[clusterID] *= real_t( n_assigned[clusterID][ALL] - 1 ); // undo previous scaling
            clusterCenter[clusterID] += barycenter[next];                         // add new point
            clusterCenter[clusterID] /= real_t( n_assigned[clusterID][ALL] );     // apply scaling
         }
      }
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
