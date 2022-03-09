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
      NIL = -1,
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
   id0[ALL + 1] = id0[ALL];                 // n_all

   // distribute number of primitives for each primitive type (round robin)
   std::vector< std::array< uint_t, ALL > > primitives_on_rank( n_processes, std::array< uint_t, ALL >{} );
   for ( auto pType = CELL; pType != NIL; pType = PT( pType - 1 ) )
   {
      for ( uint_t i = id0[pType]; i < id0[pType + 1]; ++i )
      {
         ++primitives_on_rank[i % n_processes][pType];
      }
   }

   // compute barycenter of all primitives
   std::vector< Point3D > barycenter( id0[ALL] );
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

   // which primitives are assigned to some process
   std::vector< bool > isAssigned( id0[ALL], false );
   // how many primitives of each type are assigned to each process
   std::vector< std::array< uint_t, ALL + 1 > > n_assigned( n_processes, std::array< uint_t, ALL + 1 >{} );
   // barycenter of the primitive cluster corresponding to each process
   std::vector< Point3D > clusterCenter( n_processes, Point3D() );

   // add elements
   bool done = false;
   while ( !done )
   {
      done = true;
      for ( uint_t clusterID = 0; clusterID < n_processes; ++clusterID )
      {
         // select type of next primitive

         auto pType = CELL;
         while ( pType != NIL && n_assigned[clusterID][pType] >= primitives_on_rank[clusterID][pType] )
         {
            pType = PT( pType - 1 );
         }

         if ( pType == NIL ) // no elements left to insert
         {
            continue;
         }
         else
         {
            done = false;
         }

         // range of PrimitiveIDs
         auto begin = id0[pType];
         auto end   = id0[pType + 1];
         // distribute range over processes
         auto range = ( end - begin ) / n_processes;
         auto mod   = ( end - begin ) % n_processes;
         begin += ( rank < mod ) ? ( range + 1 ) * rank : range * rank + mod;
         end = ( rank < mod ) ? begin + range + 1 : begin + range;

         auto id = begin;

         // select next primitive to add

         if ( n_assigned[clusterID][ALL] == 0 )
         {
            /* the initial element of each cluster is chosen s.th.
             id = arg max_i min_j ||barycenter[i] - clusterCenter[j]||
            */
            real_t d_max = -1.0;

            for ( uint_t i = begin; i < end; ++i )
            {
               if ( !isAssigned[i] )
               {
                  // find d_min = min_j ||barycenter[i] - clusterCenter[j]||
                  real_t d_min = std::numeric_limits< real_t >::max();
                  for ( uint_t j = 0; j < clusterID; ++j )
                  {
                     auto d = ( clusterCenter[j] - barycenter[i] ).normSq();
                     if ( d < d_min )
                     {
                        d_min = d;
                     }
                  }

                  if ( d_min > d_max )
                  {
                     d_max = d_min;
                     id    = i;
                  }
               }
            }
            auto global_max = walberla::mpi::allReduce( d_max, walberla::mpi::MAX );
            if ( global_max > d_max )
            {
               id = id0[ALL];
            }
         }
         else
         {
            // choose element with its barycenter closest to center of cluster
            real_t d_min = std::numeric_limits< real_t >::infinity();

            for ( uint_t i = begin; i < end; ++i )
            {
               if ( !isAssigned[i] )
               {
                  auto d = ( clusterCenter[clusterID] - barycenter[i] ).normSq();
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
               id = id0[ALL];
            }
         }

         id = walberla::mpi::allReduce( id, walberla::mpi::MIN );

         WALBERLA_ASSERT( id > id0[ALL] );

         // assign primitive to cluster

         if ( pType == VTX )
         {
            vtxs[id].setTargetRank( clusterID );
         }
         else if ( pType == EDGE )
         {
            edges[id - id0[EDGE]].setTargetRank( clusterID );
         }
         else if ( pType == FACE )
         {
            faces[id - id0[FACE]].setTargetRank( clusterID );
         }
         else if ( pType == CELL )
         {
            cells[id - id0[CELL]].setTargetRank( clusterID );
         }
         // mark as assigned
         ++n_assigned[clusterID][pType];
         ++n_assigned[clusterID][ALL];
         isAssigned[id] = true;
         // update center
         clusterCenter[clusterID] *= real_t( n_assigned[clusterID][ALL] - 1 ); // undo previous scaling
         clusterCenter[clusterID] += barycenter[id];                           // add new point
         clusterCenter[clusterID] /= real_t( n_assigned[clusterID][ALL] );     // apply scaling
      }
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
