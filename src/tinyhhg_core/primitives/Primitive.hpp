
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "core/NonCopyable.h"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/types/flags.hpp"

#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "core/DataTypes.h"

#include "tinyhhg_core/geometry/IdentityMap.hpp"

#include <map>
#include <memory>
#include <vector>
#include <map>

namespace hhg {

// to removed when moving to walberla namespace
using walberla::NonCopyable;
using walberla::uint_t;

namespace internal {

class PrimitiveData : private NonCopyable
{
public:

  template< typename DataType >
  PrimitiveData( const std::shared_ptr< DataType > & ptr ) : ptr_( ptr )
  {
    WALBERLA_ASSERT_NOT_NULLPTR( ptr.get() );
  }

  template< typename DataType >
  DataType* get()
  {
    return static_cast< DataType* >( ptr_.get() );
  }

private:

  std::shared_ptr< void > ptr_;

};

}

class PrimitiveID;
class PrimitiveStorage;

class Vertex;
class Edge;
class Face;
class Cell;

/// \brief Base class for primitive geometries
/// \author Nils Kohl (nils.kohl@fau.de)
///
/// The \ref Primitive class is intended to be used as a base class for primitives like vertices or edges.
///
/// Every \ref Primitive of the domain carries an unique ID a.k.a \ref PrimitiveID.
///
/// It contains methods to retrieve the \ref PrimitiveIDs from the neighboring primitives and methods
/// to (de-)serialize its metadata to / from MPI buffers.
///
/// It is able to store arbitrary data structures (e.g. from the standard library or custom classes)
/// that can, however only be added through a governing structure, for example the \ref PrimitiveStorage class.
/// Using the respective \ref PrimitiveDataID a pointer to the data can be obtained.
///
/// For more details on the data handling refer to \ref PrimitiveDataHandling and \ref PrimitiveStorage.
///
class Primitive
{
public:

  friend class SetupPrimitiveStorage;
  friend class PrimitiveStorage;

  typedef internal::PrimitiveData PrimitiveData;

  virtual ~Primitive() {}

  /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
  /// \param index the \ref PrimitiveDataID of the data that should be returned
  template< typename DataType >
  inline DataType* getData( const PrimitiveDataID< DataType, Primitive > & index ) const;

  /// Returns the number of registered data / data handling pairs.
  uint_t getNumberOfDataEntries() const { return data_.size(); }

  /// Returns the \ref PrimitiveID of the \ref Primitive
  const PrimitiveID & getID() const { return primitiveID_; }

  /// @name Neighborhood
  /// Access to IDs of neighbors of either lower or higher dimension.
  ///@{
  bool neighborPrimitiveExists( const PrimitiveID & primitiveID ) const;

  void getNeighborPrimitives( std::vector< PrimitiveID > & neighborPrimitives ) const;
  void getNeighborVertices( std::vector< PrimitiveID > & neighborVertices ) const { neighborVertices.assign( neighborVertices_.begin(), neighborVertices_.end() ); }
  void getNeighborEdges   ( std::vector< PrimitiveID > & neighborEdges )    const { neighborEdges.assign   ( neighborEdges_.begin(),    neighborEdges_.end()    ); }
  void getNeighborFaces   ( std::vector< PrimitiveID > & neighborFaces )    const { neighborFaces.assign   ( neighborFaces_.begin(),    neighborFaces_.end()    ); }
  void getNeighborCells   ( std::vector< PrimitiveID > & neighborCells )    const { neighborCells.assign   ( neighborCells_.begin(),    neighborCells_.end()    ); }

  template< typename PrimitiveType >
  inline void getNeighborPrimitivesGenerically( std::vector< PrimitiveID > & neighborPrimitives ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  const std::vector< PrimitiveID > & neighborVertices() const { return neighborVertices_; }
  const std::vector< PrimitiveID > & neighborEdges()    const { return neighborEdges_; }
  const std::vector< PrimitiveID > & neighborFaces()    const { return neighborFaces_; }
  const std::vector< PrimitiveID > & neighborCells()    const { return neighborCells_; }

  uint_t getNumNeighborPrimitives() const { return getNumNeighborVertices() + getNumNeighborEdges() + getNumNeighborFaces() + getNumNeighborCells(); }
  uint_t getNumNeighborVertices  () const { return neighborVertices_.size(); }
  uint_t getNumNeighborEdges     () const { return neighborEdges_.size(); }
  uint_t getNumNeighborFaces     () const { return neighborFaces_.size(); }
  uint_t getNumNeighborCells     () const { return neighborCells_.size(); }

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const = 0;
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const = 0;

  virtual const std::vector< PrimitiveID > & getLowerDimNeighbors()  const = 0;
  virtual const std::vector< PrimitiveID > & getHigherDimNeighbors() const = 0;

  virtual uint_t getNumLowerDimNeighbors()  const = 0;
  virtual uint_t getNumHigherDimNeighbors() const = 0;
  /// @}

  /// @name Serialization
  /// (De-)Serialize the \ref Primitive to / from MPI buffers.
  /// Even if called via the \ref Primitive class, the metadata of the
  /// subclasses (e.g. \ref Vertex or \ref Edge) are also (de-)serialized.
  /// Does not (de-)serialize attached data.
  ///@{
  void   serialize( walberla::mpi::SendBuffer & sendBuffer ) const;
  void deserialize( walberla::mpi::RecvBuffer & recvBuffer );
  ///@}

  const std::shared_ptr<GeometryMap>& getGeometryMap() const;

  /// Flag that indicates the boundary "set" that the primitive was assigned to during setup.
  /// The flag might for example originate directly from the \ref MeshInfo or can be set in the \ref SetupPrimitiveStorage.
  /// However this flag must not be used to test if a \ref Primitive matches a certain boundary condition.
  /// For this purpose use the \ref BoundaryCondition class.
  uint_t getMeshBoundaryFlag() const { return meshBoundaryFlag_; }

protected:

  /// Only subclasses shall be constructable
  /// Explicit copy-constructor - added data shall not be copied
  Primitive( const Primitive & other ) :
    primitiveID_( other.primitiveID_ ), neighborVertices_( other.neighborVertices_ ),
    neighborEdges_( other.neighborEdges_ ), neighborFaces_( other.neighborFaces_ ),
    neighborCells_( other.neighborCells_ ), geometryMap_ ( other.geometryMap_ ),
    meshBoundaryFlag_( other.meshBoundaryFlag_ )
  {}

  /// Only subclasses shall be constructable
  Primitive( const PrimitiveID & id ) :
    primitiveID_( id ), meshBoundaryFlag_( 0 ), geometryMap_( std::make_shared<IdentityMap>() )
  {}

  /// Creates Primitive from an MPI buffer
  Primitive( walberla::mpi::RecvBuffer & recvBuffer ) { deserializePrimitive( recvBuffer ); }

  template< typename DataType, typename PrimitiveType >
  inline DataType* genericGetData( const PrimitiveDataID< DataType, PrimitiveType > & index ) const;

  std::vector< PrimitiveID > neighborVertices_;
  std::vector< PrimitiveID > neighborEdges_;
  std::vector< PrimitiveID > neighborFaces_;
  std::vector< PrimitiveID > neighborCells_;

  /// Compare two primitives by their ID
  inline bool operator==(const Primitive& rhs) const {
    return this->getID() == rhs.getID();
  }

  virtual void   serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const = 0;
  virtual void deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )       = 0;

  std::shared_ptr<GeometryMap> geometryMap_;

private:

  void   serializePrimitive( walberla::mpi::SendBuffer & sendBuffer ) const;
  void deserializePrimitive( walberla::mpi::RecvBuffer & recvBuffer );

  /// Holds pointers to the attached primitive data.
  /// Mapping from the dataID index to a pointer to the data.
  std::map< uint_t, std::shared_ptr< PrimitiveData > > data_;

  PrimitiveID primitiveID_;
  uint_t      meshBoundaryFlag_;

};

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Primitive >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborPrimitives( neighborPrimitives ); }

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Vertex >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborVertices( neighborPrimitives ); }

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Edge >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborEdges( neighborPrimitives ); }

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Face >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborFaces( neighborPrimitives ); }

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Cell >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborCells( neighborPrimitives ); }


// General methods for data and data handling retrieval
template< typename DataType, typename PrimitiveType >
DataType* Primitive::genericGetData( const PrimitiveDataID< DataType, PrimitiveType > & index ) const
{
  WALBERLA_ASSERT_EQUAL( data_.count( index ), 1, "There is no data available for the specified DataID!" );
  return data_.at( index )->template get< DataType >();
}

// Methods to retrieve data and data handling from primitives
template< typename DataType >
DataType* Primitive::getData( const PrimitiveDataID< DataType, Primitive > & index ) const
{
  return genericGetData< DataType >( index );
}


} // namespace hhg

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G >   // Growth policy of SendBuffer
GenericSendBuffer<T,G>& operator<<( GenericSendBuffer<T,G> & buf, const hhg::Primitive & primitive )
{
  primitive.serialize( buf );
  return buf;
}

template< typename T >   // Element type  of RecvBuffer
GenericRecvBuffer<T>& operator>>( GenericRecvBuffer<T> & buf, hhg::Primitive & primitive )
{
  primitive.deserialize( buf );
  return buf;
}

}
}


