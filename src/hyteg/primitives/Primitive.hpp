/*
 * Copyright (c) 2017-2025 Daniel Drzisga, Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "core/DataTypes.h"
#include "core/NonCopyable.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/primitivedata/PrimitiveDataHandling.hpp"
#include "hyteg/primitivedata/PrimitiveDataID.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"
#include "hyteg/types/Concepts.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

// to be removed when moving to walberla namespace
using walberla::NonCopyable;
using walberla::uint_t;

namespace internal {

class PrimitiveData : private NonCopyable
{
 public:
   template < typename DataType >
   PrimitiveData( const std::shared_ptr< DataType >& ptr )
   : ptr_( ptr )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( ptr.get() );
   }

   template < typename DataType >
   DataType* get()
   {
      return static_cast< DataType* >( ptr_.get() );
   }

 private:
   std::shared_ptr< void > ptr_;
};

} // namespace internal

class PrimitiveID;
class PrimitiveStorage;
class SubStorage;
namespace adaptiveRefinement {
template < class K_Simplex >
class K_Mesh;
}

class Vertex;
class Edge;
class Face;
class Cell;

/// \brief Base class for primitive geometries
/// \author Nils Kohl (nils.kohl@fau.de)
///
/// The \ref Primitive class is intended to be used as a base class for primitives like vertices or edges.
///
/// Every \ref Primitive of the domain carries an unique ID a.k.a PrimitiveID.
///
/// It contains methods to retrieve the PrimitiveID from the neighboring primitives and methods
/// to (de-)serialize its metadata to / from MPI buffers.
///
/// It is able to store arbitrary data structures (e.g. from the standard library or custom classes)
/// that can, however only be added through a governing structure, for example the \ref PrimitiveStorage class.
/// Using the respective PrimitiveDataID a pointer to the data can be obtained.
///
/// For more details on the data handling refer to \ref PrimitiveDataHandling and \ref PrimitiveStorage.
///
class Primitive
{
 public:
   friend class SetupPrimitiveStorage;
   friend class PrimitiveStorage;
   friend class SubStorage;
   template < class K_Simplex >
   friend class adaptiveRefinement::K_Mesh;

   typedef internal::PrimitiveData PrimitiveData;

   virtual ~Primitive() {}

   /// Enumeration to differentiate different kinds of primitives
   enum PrimitiveTypeEnum : uint8_t
   {
      VERTEX,
      EDGE,
      FACE,
      CELL,
      INVALID
   };

   /// Return the type of the actual primitive
   virtual PrimitiveTypeEnum getType() const = 0;

   /// Returns true if the data that belongs to the passed \ref PrimitiveDataID is allocated.
   /// \param index the \ref PrimitiveDataID of the data that shall be asked for
   template < typename DataType >
   inline bool hasData( const PrimitiveDataID< DataType, Primitive >& index ) const;

   /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
   /// \param index the \ref PrimitiveDataID of the data that should be returned
   template < typename DataType >
   inline DataType* getData( const PrimitiveDataID< DataType, Primitive >& index ) const;

   /// Returns the number of registered data / data handling pairs.
   uint_t getNumberOfDataEntries() const { return data_.size(); }

   /// Returns the \ref PrimitiveID of the \ref Primitive
   const PrimitiveID& getID() const { return primitiveID_; }

   /// @name Neighborhood
   /// Access to IDs of neighbors of either lower or higher dimension.
   ///@{
   bool neighborPrimitiveExists( const PrimitiveID& primitiveID ) const;

   void getNeighborPrimitives( std::vector< PrimitiveID >& neighborPrimitives ) const;
   void getNeighborVertices( std::vector< PrimitiveID >& neighborVertices ) const
   {
      neighborVertices.assign( neighborVertices_.begin(), neighborVertices_.end() );
   }
   void getNeighborEdges( std::vector< PrimitiveID >& neighborEdges ) const
   {
      neighborEdges.assign( neighborEdges_.begin(), neighborEdges_.end() );
   }
   void getNeighborFaces( std::vector< PrimitiveID >& neighborFaces ) const
   {
      neighborFaces.assign( neighborFaces_.begin(), neighborFaces_.end() );
   }
   void getNeighborCells( std::vector< PrimitiveID >& neighborCells ) const
   {
      neighborCells.assign( neighborCells_.begin(), neighborCells_.end() );
   }

   template < typename PrimitiveType >
   inline void getNeighborPrimitivesGenerically( std::vector< PrimitiveID >& neighborPrimitives ) const
   {
      static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" );
   }

   const std::vector< PrimitiveID >& neighborVertices() const { return neighborVertices_; }
   const std::vector< PrimitiveID >& neighborEdges() const { return neighborEdges_; }
   const std::vector< PrimitiveID >& neighborFaces() const { return neighborFaces_; }
   const std::vector< PrimitiveID >& neighborCells() const { return neighborCells_; }

   uint_t getNumNeighborPrimitives() const
   {
      return getNumNeighborVertices() + getNumNeighborEdges() + getNumNeighborFaces() + getNumNeighborCells();
   }
   uint_t getNumNeighborVertices() const { return neighborVertices_.size(); }
   uint_t getNumNeighborEdges() const { return neighborEdges_.size(); }
   uint_t getNumNeighborFaces() const { return neighborFaces_.size(); }
   uint_t getNumNeighborCells() const { return neighborCells_.size(); }

   virtual void getLowerDimNeighbors( std::vector< PrimitiveID >& lowerDimNeighbors ) const   = 0;
   virtual void getHigherDimNeighbors( std::vector< PrimitiveID >& higherDimNeighbors ) const = 0;

   virtual const std::vector< PrimitiveID >& getLowerDimNeighbors() const  = 0;
   virtual const std::vector< PrimitiveID >& getHigherDimNeighbors() const = 0;

   virtual uint_t getNumLowerDimNeighbors() const  = 0;
   virtual uint_t getNumHigherDimNeighbors() const = 0;
   /// @}

   /// @name Serialization
   /// (De-)Serialize the \ref Primitive to / from MPI buffers.
   /// Even if called via the \ref Primitive class, the metadata of the
   /// subclasses (e.g. \ref Vertex or \ref Edge) are also (de-)serialized.
   /// Does not (de-)serialize attached data.
   ///@{
   void serialize( walberla::mpi::SendBuffer& sendBuffer ) const;
   void deserialize( walberla::mpi::RecvBuffer& recvBuffer );
   ///@}

   const std::shared_ptr< GeometryMap >& getGeometryMap() const;

   /// Flag that indicates the boundary "set" that the primitive was assigned to during setup.
   /// The flag might for example originate directly from the \ref MeshInfo or can be set in the \ref SetupPrimitiveStorage.
   /// However this flag must not be used to test if a \ref Primitive matches a certain boundary condition.
   /// For this purpose use the \ref BoundaryCondition class.
   uint_t getMeshBoundaryFlag() const { return meshBoundaryFlag_; }

   /// @name Refinement
   ///@{
   [[nodiscard]] bool hasChildren() const;

   [[nodiscard]] std::vector< PrimitiveID > childVertices() const;
   [[nodiscard]] std::vector< PrimitiveID > childEdges() const;
   [[nodiscard]] std::vector< PrimitiveID > childFaces() const;
   [[nodiscard]] std::vector< PrimitiveID > childCells() const;

   void addChildVertices( const std::vector< PrimitiveID >& pids );
   void addChildEdges( const std::vector< PrimitiveID >& pids );
   void addChildFaces( const std::vector< PrimitiveID >& pids );
   void addChildCells( const std::vector< PrimitiveID >& pids );

   void removeChildren( const std::vector< PrimitiveID >& pids );

   [[nodiscard]] bool hasParent() const;

   [[nodiscard]] PrimitiveID parent() const;

   void setParent( const PrimitiveID& pid );

   void clearParent();
   ///@}

 protected:
   /// Only subclasses shall be constructable
   /// Explicit copy-constructor - added data shall not be copied
   Primitive( const Primitive& other )
   : neighborVertices_( other.neighborVertices_ )
   , neighborEdges_( other.neighborEdges_ )
   , neighborFaces_( other.neighborFaces_ )
   , neighborCells_( other.neighborCells_ )
   , childVertices_( other.childVertices_ )
   , childEdges_( other.childEdges_ )
   , childFaces_( other.childFaces_ )
   , childCells_( other.childCells_ )
   , geometryMap_( other.geometryMap_ )
   , primitiveID_( other.primitiveID_ )
   , meshBoundaryFlag_( other.meshBoundaryFlag_ )
   {}

   /// Only subclasses shall be constructable
   Primitive( const PrimitiveID& id )
   : geometryMap_( std::make_shared< IdentityMap >() )
   , primitiveID_( id )
   , meshBoundaryFlag_( 0 )
   {}

   /// Creates Primitive from an MPI buffer
   Primitive( walberla::mpi::RecvBuffer& recvBuffer ) { deserializePrimitive( recvBuffer ); }

   template < typename DataType, typename PrimitiveType >
   inline bool genericHasData( const PrimitiveDataID< DataType, PrimitiveType >& index ) const;

   template < typename DataType, typename PrimitiveType >
   inline DataType* genericGetData( const PrimitiveDataID< DataType, PrimitiveType >& index ) const;

   template < typename DataType, typename PrimitiveType >
   inline void genericDeleteData( const PrimitiveDataID< DataType, PrimitiveType >& index );

   template < typename DataType >
   inline void deleteData( const PrimitiveDataID< DataType, Primitive >& index );

   std::vector< PrimitiveID > neighborVertices_;
   std::vector< PrimitiveID > neighborEdges_;
   std::vector< PrimitiveID > neighborFaces_;
   std::vector< PrimitiveID > neighborCells_;

   std::vector< PrimitiveID > childVertices_;
   std::vector< PrimitiveID > childEdges_;
   std::vector< PrimitiveID > childFaces_;
   std::vector< PrimitiveID > childCells_;

   std::optional< PrimitiveID > parent_;

   /// Compare two primitives by their ID
   inline bool operator==( const Primitive& rhs ) const { return this->getID() == rhs.getID(); }

   virtual void serializeSubclass( walberla::mpi::SendBuffer& sendBuffer ) const = 0;
   virtual void deserializeSubclass( walberla::mpi::RecvBuffer& recvBuffer )     = 0;

   std::shared_ptr< GeometryMap > geometryMap_;

 private:
   void serializePrimitive( walberla::mpi::SendBuffer& sendBuffer ) const;
   void deserializePrimitive( walberla::mpi::RecvBuffer& recvBuffer );

   /// Holds pointers to the attached primitive data.
   /// Mapping from the dataID index to a pointer to the data.
   std::map< uint_t, std::shared_ptr< PrimitiveData > > data_;

   PrimitiveID primitiveID_;
   uint_t      meshBoundaryFlag_;
};

template <>
inline void Primitive::getNeighborPrimitivesGenerically< Primitive >( std::vector< PrimitiveID >& neighborPrimitives ) const
{
   getNeighborPrimitives( neighborPrimitives );
}

template <>
inline void Primitive::getNeighborPrimitivesGenerically< Vertex >( std::vector< PrimitiveID >& neighborPrimitives ) const
{
   getNeighborVertices( neighborPrimitives );
}

template <>
inline void Primitive::getNeighborPrimitivesGenerically< Edge >( std::vector< PrimitiveID >& neighborPrimitives ) const
{
   getNeighborEdges( neighborPrimitives );
}

template <>
inline void Primitive::getNeighborPrimitivesGenerically< Face >( std::vector< PrimitiveID >& neighborPrimitives ) const
{
   getNeighborFaces( neighborPrimitives );
}

template <>
inline void Primitive::getNeighborPrimitivesGenerically< Cell >( std::vector< PrimitiveID >& neighborPrimitives ) const
{
   getNeighborCells( neighborPrimitives );
}

// General methods for data and data handling retrieval
template < typename DataType, typename PrimitiveType >
bool Primitive::genericHasData( const PrimitiveDataID< DataType, PrimitiveType >& index ) const
{
   return data_.count( index ) > 0;
}

template < typename DataType, typename PrimitiveType >
DataType* Primitive::genericGetData( const PrimitiveDataID< DataType, PrimitiveType >& index ) const
{
   WALBERLA_ASSERT( genericHasData( index ), "Cannot retrieve data - there is no data available for the specified DataID!" );
   return data_.at( index )->template get< DataType >();
}

template < typename DataType, typename PrimitiveType >
void Primitive::genericDeleteData( const PrimitiveDataID< DataType, PrimitiveType >& index )
{
   data_.erase( index );
}

// Methods to retrieve data and data handling from primitives
template < typename DataType >
bool Primitive::hasData( const PrimitiveDataID< DataType, Primitive >& index ) const
{
   return genericHasData< DataType >( index );
}

template < typename DataType >
DataType* Primitive::getData( const PrimitiveDataID< DataType, Primitive >& index ) const
{
   return genericGetData< DataType >( index );
}

template < typename DataType >
void Primitive::deleteData( const PrimitiveDataID< DataType, Primitive >& index )
{
   return genericDeleteData< DataType >( index );
}

} // namespace hyteg

namespace walberla {
namespace mpi {

template < typename T,  // Element type of SendBuffer
           typename G > // Growth policy of SendBuffer
GenericSendBuffer< T, G >& operator<<( GenericSendBuffer< T, G >& buf, const hyteg::Primitive& primitive )
{
   primitive.serialize( buf );
   return buf;
}

template < typename T > // Element type  of RecvBuffer
GenericRecvBuffer< T >& operator>>( GenericRecvBuffer< T >& buf, hyteg::Primitive& primitive )
{
   primitive.deserialize( buf );
   return buf;
}

} // namespace mpi
} // namespace walberla
