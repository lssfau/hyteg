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

void inheritRankFromVolumePrimitives( std::map< PrimitiveID, VertexData >& vtxs,
                                      std::map< PrimitiveID, EdgeData >&   edges,
                                      std::map< PrimitiveID, FaceData >&   faces,
                                      std::map< PrimitiveID, CellData >&   cells,
                                      const NeighborhoodMap&               nbrHood )
{
   const PrimitiveType VOL = ( cells.size() == 0 ) ? FACE : CELL;

   // loop over all interface primitive types
   for ( auto pt = VTX; pt != VOL; pt = PrimitiveType( pt + 1 ) )
   {
      // loop over all interfaces of type pt
      for ( auto& [interfaceID, interfaceNbrs] : nbrHood[pt] )
      {
         // loop over all neighboring volume primitives and extract their target ranks
         std::map< uint_t, uint_t > nbrRnks;
         for ( auto& nbrID : interfaceNbrs[VOL] )
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
   }
}

} // namespace adaptiveRefinement
} // namespace hyteg
