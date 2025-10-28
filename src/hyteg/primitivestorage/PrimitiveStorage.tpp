/*
 * Copyright (c) 2017-2025 Dominik Bartuschat, Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

template <>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Primitive >( const PrimitiveID& id ) const
{
   return primitiveExistsLocally( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Vertex >( const PrimitiveID& id ) const
{
   return vertexExistsLocally( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Edge >( const PrimitiveID& id ) const
{
   return edgeExistsLocally( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Face >( const PrimitiveID& id ) const
{
   return faceExistsLocally( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsLocallyGenerically< Cell >( const PrimitiveID& id ) const
{
   return cellExistsLocally( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Primitive >( const PrimitiveID& id ) const
{
   return primitiveExistsInNeighborhood( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Vertex >( const PrimitiveID& id ) const
{
   return vertexExistsInNeighborhood( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Edge >( const PrimitiveID& id ) const
{
   return edgeExistsInNeighborhood( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Face >( const PrimitiveID& id ) const
{
   return faceExistsInNeighborhood( id );
}

template <>
inline bool PrimitiveStorage::primitiveExistsInNeighborhoodGenerically< Cell >( const PrimitiveID& id ) const
{
   return cellExistsInNeighborhood( id );
}

template <>
inline const Primitive* PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID& id ) const
{
   return getPrimitive( id );
}

template <>
inline Primitive* PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID& id )
{
   return getPrimitive( id );
}

template <>
inline const Vertex* PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID& id ) const
{
   return getVertex( id );
}

template <>
inline Vertex* PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID& id )
{
   return getVertex( id );
}

template <>
inline const Edge* PrimitiveStorage::getPrimitiveGenerically< Edge >( const PrimitiveID& id ) const
{
   return getEdge( id );
}

template <>
inline Edge* PrimitiveStorage::getPrimitiveGenerically< Edge >( const PrimitiveID& id )
{
   return getEdge( id );
}

template <>
inline const Face* PrimitiveStorage::getPrimitiveGenerically< Face >( const PrimitiveID& id ) const
{
   return getFace( id );
}

template <>
inline Face* PrimitiveStorage::getPrimitiveGenerically< Face >( const PrimitiveID& id )
{
   return getFace( id );
}

template <>
inline const Cell* PrimitiveStorage::getPrimitiveGenerically< Cell >( const PrimitiveID& id ) const
{
   return getCell( id );
}

template <>
inline Cell* PrimitiveStorage::getPrimitiveGenerically< Cell >( const PrimitiveID& id )
{
   return getCell( id );
}

template <>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Primitive >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getPrimitiveIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Vertex >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getVertexIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Edge >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getEdgeIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Face >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getFaceIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getPrimitiveIDsGenerically< Cell >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getCellIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getNeighboringPrimitiveIDsGenerically< Vertex >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getNeighboringVertexIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getNeighboringPrimitiveIDsGenerically< Edge >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getNeighboringEdgeIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getNeighboringPrimitiveIDsGenerically< Face >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getNeighboringFaceIDs( primitiveIDs );
}

template <>
inline void PrimitiveStorage::getNeighboringPrimitiveIDsGenerically< Cell >( std::vector< PrimitiveID >& primitiveIDs ) const
{
   getNeighboringCellIDs( primitiveIDs );
}
