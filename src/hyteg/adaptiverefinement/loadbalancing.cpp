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
      VTX,
      EDGE,
      FACE,
      CELL,
      ALL
   };
   /* we assume that the elements in the input vectors are ordered by PrimitiveID
      and that for each vertex v, edge e, face f and cell c it holds
               id_v < id_e < id_f < id_c
   */
   std::array< uint_t, 6 > id0;
   id0[VTX]  = 0;                        // vtxID0
   id0[EDGE] = id0[VTX] + vtxs.size();   // edgeID0
   id0[FACE] = id0[EDGE] + edges.size(); // faceID0
   id0[CELL] = id0[FACE] + faces.size(); // cellID0
   id0[ALL]  = id0[CELL] + cells.size(); // n_all
   id0[5]    = id0[ALL];                 // n_all

   // distribute number of primitives for each primitive type (round robin)
   std::vector< std::array< uint_t, 4 > > primitives_on_rank( n_processes, std::array< uint_t, 4 >{} );
   for ( auto pType = VTX; pType < ALL; pType = PT( pType + 1 ) )
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

   // find an appropriate cluster of primitives for each rank
   std::vector< bool > isAssigned( id0[ALL], false );
   for ( uint_t clusterID = 0; clusterID < n_processes; ++clusterID )
   {
      // number of primitives assigned to this cluster [v, e, f, c, all]
      std::array< uint_t, 5 > n_assigned{};
      // center of cluster
      Point3D center;

      // choose first element
      uint_t id = 0;
      while ( id < id0[ALL] && isAssigned[id] )
      {
         ++id;
      }
      // type of primitive
      auto pType = CELL;
      while ( id < id0[pType] )
      {
         pType = PT( pType - 1 );
      }

      // add elements
      while ( id < id0[ALL] )
      {
         // assign primitive to this cluster
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
         ++n_assigned[pType];
         ++n_assigned[ALL];
         isAssigned[id] = true;

         // update center
         center *= real_t( n_assigned[ALL] - 1 ); // undo previous scaling
         center += barycenter[id];                // add new point
         center /= real_t( n_assigned[ALL] );     // apply scaling

         // choose type of next element
         pType = VTX;
         while ( pType < ALL && n_assigned[pType] >= primitives_on_rank[clusterID][pType] )
         {
            pType = PT( pType + 1 );
         }
         auto begin = id0[pType];
         auto end   = id0[pType + 1];

         // parallelize
         auto range = ( end - begin ) / n_processes;
         auto mod   = ( end - begin ) % n_processes;
         begin += ( rank < mod ) ? ( range + 1 ) * rank : range * rank + mod;
         end = ( rank < mod ) ? begin + range + 1 : begin + range;

         // choose element with its barycenter closest to center of cluster (MPI parallel)
         real_t d_min = std::numeric_limits< real_t >::infinity();
         id           = id0[ALL];
         for ( uint_t i = begin; i < end; ++i )
         {
            if ( !isAssigned[i] )
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
            id = id0[ALL];
         }
         id = walberla::mpi::allReduce( id, walberla::mpi::MIN );
      }
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
