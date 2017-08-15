
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "core/NonCopyable.h"
#include "tinyhhg_core/primitiveid.hpp"

#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"

#include <memory>
#include <vector>

namespace hhg {

// to removed when moving to walberla namespace
using walberla::NonCopyable;

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
#if 0
class ConstPrimitiveData : private NonCopyable
{
public:

  template< typename DataType >
  ConstPrimitiveData( const std::shared_ptr< const DataType > & ptr ) : ptr_( ptr )
  {
    WALBERLA_ASSERT_NOT_NULLPTR( ptr.get() );
  }

  template< typename DataType >
  const DataType* get()
  {
    return static_cast< const DataType* >( ptr_.get() );
  }

private:

  std::shared_ptr< const void > ptr_;

};
#endif
}

class PrimitiveID;
class PrimitiveStorage;

class Vertex;
class Edge;
class Face;

/// \brief Base class for primitive geometries
/// \author Nils Kohl (nils.kohl@fau.de)
///
/// The \ref Primitive class is intended to be used as a base class for primitives like vertices or edges.
///
/// It is able to store arbitrary data structures (e.g. from the standard library or custom classes)
/// that can, however only be added through a governing structure, for example the \ref PrimitiveStorage class.
/// Using the respective \ref PrimitiveDataID a pointer to the data or the associated \ref PrimitiveDataHandling can be obtained.
///
/// For more details on the data handling refer to \ref PrimitiveDataHandling.
///
class Primitive
{
public:

  friend class PrimitiveStorage;

  typedef internal::PrimitiveData PrimitiveData;
#if 0
  typedef internal::ConstPrimitiveData ConstPrimitiveData;
#endif

  virtual ~Primitive() {}

  /// Returns a pointer to the data that belongs to the passed \ref PrimitiveDataID.
  /// \param index the \ref PrimitiveDataID of the data that should be returned
  template< typename DataType >
  inline DataType* getData( const PrimitiveDataID< DataType, Primitive > & index ) const;

  /// Returns a pointer to the \ref PrimitiveDataHandling that belongs to the passed \ref PrimitiveDataID.
  /// \param index the \ref PrimitiveDataID of the data handling that should be returned
  template< typename DataType >
  inline PrimitiveDataHandling< DataType, Primitive >* getDataHandling( const PrimitiveDataID< DataType, Primitive > & index ) const;

  /// Returns the number of registered data / data handling pairs.
  uint_t getNumberOfDataEntries() const { return data_.size(); }

  const PrimitiveID & getID() const { return primitiveID_; }

  /// @name Neighborhood
  /// Access to IDs of neighbors of either lower or higher dimension.
  ///@{
  void getNeighborPrimitives( std::vector< PrimitiveID > & neighborPrimitives ) const;
  void getNeighborVertices( std::vector< PrimitiveID > & neighborVertices ) const { neighborVertices.assign( neighborVertices_.begin(), neighborVertices_.end() ); }
  void getNeighborEdges   ( std::vector< PrimitiveID > & neighborEdges )    const { neighborEdges.assign   ( neighborEdges_.begin(),    neighborEdges_.end()    ); }
  void getNeighborFaces   ( std::vector< PrimitiveID > & neighborFaces )    const { neighborFaces.assign   ( neighborFaces_.begin(),    neighborFaces_.end()    ); }

  template< typename PrimitiveType >
  inline void getNeighborPrimitivesGenerically( std::vector< PrimitiveID > & neighborPrimitives ) const { static_assert( sizeof( PrimitiveType ) == 0 /* always false */, "Invalid primitive type" ); }

  const std::vector< PrimitiveID > & neighborVertices() const { return neighborVertices_; }
  const std::vector< PrimitiveID > & neighborEdges()    const { return neighborEdges_; }
  const std::vector< PrimitiveID > & neighborFaces()    const { return neighborFaces_; }

  uint_t getNumNeighborVertices() const { return neighborVertices_.size(); }
  uint_t getNumNeighborEdges   () const { return neighborEdges_.size(); }
  uint_t getNumNeighborFaces   () const { return neighborFaces_.size(); }

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const = 0;
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const = 0;

  virtual const std::vector< PrimitiveID > & getLowerDimNeighbors()  const = 0;
  virtual const std::vector< PrimitiveID > & getHigherDimNeighbors() const = 0;

  virtual uint_t getNumLowerDimNeighbors()  const = 0;
  virtual uint_t getNumHigherDimNeighbors() const = 0;
  /// @}

protected:

  /// Only subclasses shall be constructable
  /// Explicit copy-constructor - added data shall not be copied
  Primitive( const Primitive & other ) :
    primitiveID_( other.primitiveID_ ), neighborVertices_( other.neighborVertices_ ),
    neighborEdges_( other.neighborEdges_ ), neighborFaces_( other.neighborFaces_ )
  {}

  /// Only subclasses shall be constructable
  Primitive( const PrimitiveID & id ) : primitiveID_( id ) {}

  template< typename DataType, typename PrimitiveType >
  inline DataType* genericGetData( const PrimitiveDataID< DataType, PrimitiveType > & index ) const;

#if 0
  template< typename DataType, typename PrimitiveType >
  inline PrimitiveDataHandling< DataType, PrimitiveType >* getDataHandling( const PrimitiveDataID< DataType, PrimitiveType > & index ) const;

  /// Not public in order to guarantee that data is only added through the governing structure.
  /// This ensures valid DataIDs.
  template< typename DataType, typename DataHandlingType >
  inline void addData( const PrimitiveDataID< DataType, Primitive > & index,
                       const std::shared_ptr< DataHandlingType > & dataHandling );

  template< typename DataType, typename PrimitiveType, typename DataHandlingType >
  inline void genericAddData( const PrimitiveDataID< DataType, PrimitiveType > & index,
                              const std::shared_ptr< DataHandlingType > & dataHandling,
		                          const PrimitiveType * primitive );
#endif

  std::vector< PrimitiveID > neighborVertices_;
  std::vector< PrimitiveID > neighborEdges_;
  std::vector< PrimitiveID > neighborFaces_;

  /// Compare two primitives by their ID
  inline bool operator==(const Primitive& rhs) const {
    return this->getID() == rhs.getID();
  }

private:

  /// Holds pointers to the attached primitive data.
  /// Mapping from the dataID index to a pointer to the data.
  std::map< uint_t, std::shared_ptr< PrimitiveData > > data_;

  PrimitiveID primitiveID_;

};

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Primitive >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborPrimitives( neighborPrimitives ); }

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Vertex >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborVertices( neighborPrimitives ); }

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Edge >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborEdges( neighborPrimitives ); }

template<>
inline void Primitive::getNeighborPrimitivesGenerically< Face >( std::vector< PrimitiveID > & neighborPrimitives ) const { getNeighborFaces( neighborPrimitives ); }

// General methods for data and data handling retrieval
template< typename DataType, typename PrimitiveType >
DataType* Primitive::genericGetData( const PrimitiveDataID< DataType, PrimitiveType > & index ) const
{
  WALBERLA_ASSERT_EQUAL( data_.count( index ), 1, "There is no data available for the specified index" );
  return data_.at( index )->template get< DataType >();
}

#if 0
template< typename DataType, typename PrimitiveType >
PrimitiveDataHandling< DataType, PrimitiveType >* Primitive::getDataHandling( const PrimitiveDataID< DataType, PrimitiveType > & index ) const
{
  WALBERLA_ASSERT_EQUAL( data_.count( index ), 1, "There is no data handling available for the specified index" );
  return data_.at( index ).second->template get< PrimitiveDataHandling< DataType, PrimitiveType > >();
}
///////////////////////////////////////////////////////
#endif

// Methods to retrieve data and data handling from primitives
template< typename DataType >
DataType* Primitive::getData( const PrimitiveDataID< DataType, Primitive > & index ) const
{
  return genericGetData< DataType >( index );
}

#if 0
template< typename DataType >
PrimitiveDataHandling< DataType, Primitive >* Primitive::getDataHandling( const PrimitiveDataID< DataType, Primitive > & index ) const
{
  WALBERLA_ASSERT_EQUAL( data_.count( index ), 1, "There is no data handling available for the specified index" );
  return data_.at( index ).second->template get< PrimitiveDataHandling< DataType, Primitive > >();
}
/////////////////////////////////////////////////////////////
#endif

#if 0
template< typename DataType, typename DataHandlingType >
void Primitive::addData( const PrimitiveDataID< DataType, Primitive > & index,
                         const std::shared_ptr< DataHandlingType > & dataHandling )
{
  genericAddData( index, dataHandling, this );
}


template< typename DataType,
          typename PrimitiveType,
          typename DataHandlingType >
void Primitive::genericAddData( const PrimitiveDataID< DataType, PrimitiveType > & index,
                                const std::shared_ptr< DataHandlingType > & dataHandling,
                                const PrimitiveType * const primitive )
{
  WALBERLA_ASSERT_NOT_NULLPTR( dataHandling.get() );
  WALBERLA_ASSERT_NOT_NULLPTR( primitive );

  WALBERLA_ASSERT_EQUAL( data_.count( index ), 0,                        "Data at DataID index already exists" );
  WALBERLA_ASSERT_EQUAL( dataInitializationFunctions_.count( index ), 0, "Init callback at DataID index already exists" );

  WALBERLA_ASSERT_EQUAL( data_.size(), dataInitializationFunctions_.size() );

  // Set up initialization, serialization and deserialization callbacks
  auto initCallback = [ this, index, dataHandling, primitive ]() -> void
  {
    data_[index].first = std::shared_ptr< PrimitiveData >( new PrimitiveData( dataHandling->initialize( primitive ) ) );
  };

  std::function< void( walberla::mpi::SendBuffer & sendBuffer ) > serializationCallback = [ dataHandling, primitive, index ]( walberla::mpi::SendBuffer & sendBuffer ) -> void
  {
    dataHandling->serialize( primitive, index, sendBuffer );
  };

  std::function< void( walberla::mpi::RecvBuffer & recvBuffer ) > deserializationCallback = [ dataHandling, primitive, index ]( walberla::mpi::RecvBuffer & recvBuffer ) -> void
  {
    dataHandling->deserialize( primitive, index, recvBuffer );
  };

  dataInitializationFunctions_[index]  = initCallback;
  dataSerializationFunctions_[index]   = serializationCallback;
  dataDeserializationFunctions_[index] = deserializationCallback;

  initCallback();
}
#endif

} // namespace hhg

