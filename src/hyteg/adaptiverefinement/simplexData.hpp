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

#pragma once

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "simplex.hpp"

namespace hyteg {
namespace adaptiveRefinement {

enum Locality
{
   NONE  = 0,
   HALO  = 1,
   LOCAL = 2
};

// stores only essential data of a simplex to create PrimitiveStorage
template < uint_t J >
class SimplexData
{
 public:
   SimplexData( PrimitiveID                 geometryMap  = PrimitiveID(),
                uint_t                      boundaryFlag = 0,
                PrimitiveID                 id           = PrimitiveID(),
                std::array< uint_t, J + 1 > vtxIncices   = {},
                uint_t                      targetRank   = 0 )
   : _vertices( vtxIncices )
   , _geometryMap( geometryMap )
   , _boundaryFlag( boundaryFlag )
   , _id( id )
   , _targetRank( targetRank )
   , _locality( NONE )
   {}

   template < class J_Simplex >
   SimplexData( const Simplex< J, J_Simplex >* simplex )
   : _vertices( simplex->get_vertices() )
   , _geometryMap( simplex->getGeometryMap() )
   , _boundaryFlag( simplex->getBoundaryFlag() )
   , _id( simplex->getPrimitiveID() )
   , _targetRank( simplex->getTargetRank() )
   , _locality( NONE )
   {}

   inline const std::array< uint_t, J + 1 >& get_vertices() const { return _vertices; }
   inline const PrimitiveID&                 getGeometryMap() const { return _geometryMap; }
   inline const uint_t&                      getBoundaryFlag() const { return _boundaryFlag; }
   inline const PrimitiveID&                 getPrimitiveID() const { return _id; }
   inline const uint_t&                      getTargetRank() const { return _targetRank; }
   inline void                               setTargetRank( uint_t rnk ) { _targetRank = rnk; }
   inline bool                               onHalo() const { return _locality == HALO; }
   inline bool                               isLocal() const { return _locality == LOCAL; }
   inline void                               setLocality( Locality loc ) { _locality = loc; }

   // count number of vertices that are not shared between this and other, i.e.,
   // return d = (J+1) - n_shared
   // e.g., if *this==other, then d=0; if other is a direct neighbor of this, then d=1;
   inline uint_t diff( const SimplexData& other ) const
   {
      uint_t d = 0;
      for ( auto& v : _vertices )
      {
         if ( std::find( other.get_vertices().begin(), other.get_vertices().end(), v ) == other.get_vertices().end() )
            ++d;
      }
      return d;
   }

   inline void serialize( walberla::mpi::SendBuffer& sendBuffer ) const
   {
      sendBuffer << _vertices;
      sendBuffer << _geometryMap;
      sendBuffer << _boundaryFlag;
      sendBuffer << _id;
      sendBuffer << _targetRank;
   }
   inline void deserialize( walberla::mpi::RecvBuffer& recvBuffer )
   {
      recvBuffer >> _vertices;
      recvBuffer >> _geometryMap;
      recvBuffer >> _boundaryFlag;
      recvBuffer >> _id;
      recvBuffer >> _targetRank;
   }

   inline std::array< Point3D, J + 1 > get_coordinates( const EnumeratedList< Point3D >& global_coordinates ) const
   {
      std::array< Point3D, J + 1 > coords;
      for ( uint_t i = 0; i <= J; ++i )
      {
         coords[i] = global_coordinates[_vertices[i]];
      }
      return coords;
   }

 private:
   std::array< uint_t, J + 1 > _vertices;
   PrimitiveID                 _geometryMap;
   uint_t                      _boundaryFlag;
   PrimitiveID                 _id;
   uint_t                      _targetRank;
   Locality                    _locality;
};

using VertexData = SimplexData< VTX >;
using EdgeData   = SimplexData< EDGE >;
using FaceData   = SimplexData< FACE >;
using CellData   = SimplexData< CELL >;

} // namespace adaptiveRefinement
} // namespace hyteg

namespace walberla {
namespace mpi {

template < uint_t J,    // spatial dimension
           typename T,  // Element type of SendBuffer
           typename G > // Growth policy of SendBuffer
GenericSendBuffer< T, G >& operator<<( GenericSendBuffer< T, G >& buf, const hyteg::adaptiveRefinement::SimplexData< J >& sd )
{
   sd.serialize( buf );
   return buf;
}

template < uint_t J,    // spatial dimension
           typename T > // Element type  of RecvBuffer
GenericRecvBuffer< T >& operator>>( GenericRecvBuffer< T >& buf, hyteg::adaptiveRefinement::SimplexData< J >& sd )
{
   sd.deserialize( buf );
   return buf;
}

} // namespace mpi
} // namespace walberla