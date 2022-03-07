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
// #include <core/Environment.h>
// #include <core/Format.hpp>
// #include <core/logging/all.h>
#include <core/mpi/Broadcast.h>
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
   /* we assume that the elements in the input vectors are ordered by PrimitiveID
      and that for each vertex v, edge e, face f and cell c it holds
               id_v < id_e < id_f < id_c
   */
   const uint_t edgeID0 = vtxs.size();
   const uint_t faceID0 = edgeID0 + edges.size();
   const uint_t cellID0 = faceID0 + faces.size();
   const uint_t n_all   = cellID0 + cells.size();

   std::vector< uint_t > vtxs_on_rnk( n_processes, 0 );
   std::vector< uint_t > edges_on_rnk( n_processes, 0 );
   std::vector< uint_t > faces_on_rnk( n_processes, 0 );
   std::vector< uint_t > cells_on_rnk( n_processes, 0 );

   // distribute number of primitives for each primitive type
   for ( uint_t i = 0; i < edgeID0; ++i )
   {
      ++vtxs_on_rnk[i % n_processes];
   }
   for ( uint_t i = edgeID0; i < faceID0; ++i )
   {
      ++edges_on_rnk[i % n_processes];
   }
   for ( uint_t i = faceID0; i < cellID0; ++i )
   {
      ++faces_on_rnk[i % n_processes];
   }
   for ( uint_t i = cellID0; i < n_all; ++i )
   {
      ++cells_on_rnk[i % n_processes];
   }

   // compute barycenter of all primitives
   std::vector< Point3D > barycenter( n_all );
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

   // find an appropriate cluster of primitives for each rank
   std::vector< bool > hasRank( n_all, false );
   for ( uint_t clusterID = 0; clusterID < n_processes; ++clusterID )
   {
      uint_t vtxs_in_cluster  = 0;
      uint_t edges_in_cluster = 0;
      uint_t faces_in_cluster = 0;
      uint_t cells_in_cluster = 0;
      uint_t cluster_size     = 0;

      // center of cluster
      Point3D center;

      // choose first element
      uint_t id = 0;
      while ( hasRank[id] && id < n_all )
      {
         ++id;
      }
      // type of primitive
      uint_t pType = 0;
      if ( id >= edgeID0 )
      {
         pType = 1;
      }
      if ( id >= faceID0 )
      {
         pType = 2;
      }
      if ( id >= cellID0 )
      {
         pType = 3;
      }

      // add elements
      while ( id < n_all )
      {
         // add element
         if ( pType == 0 )
         {
            vtxs[id].setTargetRank( clusterID );
            ++vtxs_in_cluster;
         }
         else if ( pType == 1 )
         {
            edges[id - edgeID0].setTargetRank( clusterID );
            ++edges_in_cluster;
         }
         else if ( pType == 2 )
         {
            faces[id - faceID0].setTargetRank( clusterID );
            ++faces_in_cluster;
         }
         else if ( pType == 3 )
         {
            cells[id - cellID0].setTargetRank( clusterID );
            ++cells_in_cluster;
         }

         // mark as added
         hasRank[id] = true;

         // update center
         center *= real_t( cluster_size ); // undo previous scaling
         center += barycenter[id];         // add new point
         ++cluster_size;                   // increase cluster size
         center /= real_t( cluster_size ); // apply scaling

         // choose type of next element
         uint_t begin = n_all;
         uint_t end   = n_all;
         if ( vtxs_in_cluster < vtxs_on_rnk[clusterID] )
         {
            pType = 0;
            begin = 0;
            end   = edgeID0;
         }
         else if ( edges_in_cluster < edges_on_rnk[clusterID] )
         {
            pType = 1;
            begin = edgeID0;
            end   = faceID0;
         }
         else if ( faces_in_cluster < faces_on_rnk[clusterID] )
         {
            pType = 2;
            begin = faceID0;
            end   = cellID0;
         }
         else if ( cells_in_cluster < cells_on_rnk[clusterID] )
         {
            pType = 3;
            begin = cellID0;
            end   = n_all;
         }

         // distribute range of possible elements over processes
         auto range_tot = end - begin;
         auto range     = range_tot / n_processes;
         auto mod       = range_tot % n_processes;
         begin += ( rank < mod ) ? ( range + 1 ) * rank : range * rank + mod;
         end = ( rank < mod ) ? begin + range + 1 : begin + range;

         // choose element with its barycenter closest to center of cluster (MPI parallel)
         real_t d_min = std::numeric_limits< real_t >::infinity();
         id           = n_all;
         for ( uint_t i = begin; i < end; ++i )
         {
            if ( !hasRank[i] )
            {
               auto d = ( center - barycenter[i] ).normSq();
               if ( d < d_min )
               {
                  d_min = d;
                  id    = i;
               }
            }
         }
         auto global_min = walberla::mpi::allReduce( d_min, walberla::mpi::MIN );
         if ( global_min < d_min )
         {
            id = n_all;
         }
         id = walberla::mpi::allReduce( id, walberla::mpi::MIN );
      }
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
