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
inline std::shared_ptr< const Primitive > PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID& id ) const
{
   return getPointerToPrimitive( id );
}

template <>
inline std::shared_ptr< Primitive > PrimitiveStorage::getPrimitiveGenerically< Primitive >( const PrimitiveID& id )
{
   return getPointerToPrimitive( id );
}

template <>
inline std::shared_ptr< const Vertex > PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID& id ) const
{
   return getPointerToVertex( id );
}

template <>
inline std::shared_ptr< Vertex > PrimitiveStorage::getPrimitiveGenerically< Vertex >( const PrimitiveID& id )
{
   return getPointerToVertex( id );
}

template <>
inline std::shared_ptr< const Edge > PrimitiveStorage::getPrimitiveGenerically< Edge >( const PrimitiveID& id ) const
{
   return getPointerToEdge( id );
}

template <>
inline std::shared_ptr< Edge > PrimitiveStorage::getPrimitiveGenerically< Edge >( const PrimitiveID& id )
{
   return getPointerToEdge( id );
}

template <>
inline std::shared_ptr< const Face > PrimitiveStorage::getPrimitiveGenerically< Face >( const PrimitiveID& id ) const
{
   return getPointerToFace( id );
}

template <>
inline std::shared_ptr< Face > PrimitiveStorage::getPrimitiveGenerically< Face >( const PrimitiveID& id )
{
   return getPointerToFace( id );
}

template <>
inline std::shared_ptr< const Cell > PrimitiveStorage::getPrimitiveGenerically< Cell >( const PrimitiveID& id ) const
{
   return getPointerToCell( id );
}

template <>
inline std::shared_ptr< Cell > PrimitiveStorage::getPrimitiveGenerically< Cell >( const PrimitiveID& id )
{
   return getPointerToCell( id );
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
