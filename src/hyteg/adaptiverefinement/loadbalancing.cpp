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
MigrationInfo loadbalancing( std::map< PrimitiveID, VertexData >& vtxs,
                             std::map< PrimitiveID, EdgeData >&   edges,
                             std::map< PrimitiveID, FaceData >&   faces,
                             std::map< PrimitiveID, CellData >&   cells,
                             const uint_t&                        n_processes,
                             const uint_t&                        rank )
{
   MigrationMap_T migrationMap;
   uint_t         numReceivingPrimitives = 0;

   // roundrobin
   uint_t i = 0;

   for ( auto& [id, vtx] : vtxs )
   {
      auto currentRnk = vtx.getTargetRank();
      auto targetRnk  = i % n_processes;

      vtx.setTargetRank( targetRnk );

      if ( rank == currentRnk )
      {
         migrationMap[id] = targetRnk;
      }
      if ( rank == targetRnk )
      {
         ++numReceivingPrimitives;
      }

      ++i;
   }
   for ( auto& [id, edge] : edges )
   {
      auto currentRnk = edge.getTargetRank();
      auto targetRnk  = i % n_processes;

      edge.setTargetRank( targetRnk );

      if ( rank == currentRnk )
      {
         migrationMap[id] = targetRnk;
      }
      if ( rank == targetRnk )
      {
         ++numReceivingPrimitives;
      }

      ++i;
   }
   for ( auto& [id, face] : faces )
   {
      auto currentRnk = face.getTargetRank();
      auto targetRnk  = i % n_processes;

      face.setTargetRank( targetRnk );

      if ( rank == currentRnk )
      {
         migrationMap[id] = targetRnk;
      }
      if ( rank == targetRnk )
      {
         ++numReceivingPrimitives;
      }

      ++i;
   }
   for ( auto& [id, cell] : cells )
   {
      auto currentRnk = cell.getTargetRank();
      auto targetRnk  = i % n_processes;

      cell.setTargetRank( targetRnk );

      if ( rank == currentRnk )
      {
         migrationMap[id] = targetRnk;
      }
      if ( rank == targetRnk )
      {
         ++numReceivingPrimitives;
      }

      ++i;
   }

   return MigrationInfo( migrationMap, numReceivingPrimitives );
}

void inheritRankFromVolumePrimitives( std::map< PrimitiveID, VertexData >&         vtxs,
                                      std::map< PrimitiveID, EdgeData >&           edges,
                                      std::map< PrimitiveID, FaceData >&           faces,
                                      std::map< PrimitiveID, CellData >&           cells,
                                      const std::map< PrimitiveID, Neighborhood >& nbrHood )
{
   using PT     = PrimitiveType;
   const PT VOL = ( cells.size() == 0 ) ? FACE : CELL;

   // compute neighboring volume primitives of all interface primitives
   std::array< std::map< PrimitiveID, std::vector< PrimitiveID > >, PrimitiveType::ALL > nbrVolumes;
   for ( auto& [volID, volNbrHood] : nbrHood )
   {
      auto pt = VTX;
      while ( pt != VOL )
      {
         for ( auto& nbrID : volNbrHood[pt] )
         {
            nbrVolumes[pt][nbrID].push_back( volID );
         }

         pt = PT( pt + 1 );
      }
   }

   // loop over all interface primitive types
   auto pt = VTX;
   while ( pt != VOL )
   {
      // loop over all interfaces of type pt
      for ( auto& [interfaceID, myNbrVolumes] : nbrVolumes[pt] )
      {
         // loop over all neighboring volume primitives and extract their target ranks
         std::map< uint_t, uint_t > nbrRnks;
         for ( auto& nbrID : myNbrVolumes )
         {
            auto targetRank = ( VOL == CELL ) ? cells[nbrID].getTargetRank() : faces[nbrID].getTargetRank();
            nbrRnks[targetRank]++;
         }

         // find rank where the most neighboring volume primitives are located
         auto cmp = []( const std::pair< uint_t, uint_t >& a, const std::pair< uint_t, uint_t >& b ) {
            return a.second < b.second;
         };
         auto targetRank = std::max_element( nbrRnks.begin(), nbrRnks.end(), cmp )->first;

         // assign interface to targetRank
         switch ( pt )
         {
         case VTX:
            vtxs[interfaceID].setTargetRank( targetRank );
            break;
         case EDGE:
            edges[interfaceID].setTargetRank( targetRank );
            break;
         case FACE:
            faces[interfaceID].setTargetRank( targetRank );
            break;
         default:
            break;
         }
      }

      pt = PT( pt + 1 );
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
