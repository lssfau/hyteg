
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "core/NonCopyable.h"
#include "tinyhhg_core/primitiveid.hpp"
#include "tinyhhg_core/primitives/SetupPrimitive.hpp"

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
  PrimitiveData( DataType* ptr ) : ptr_( ptr )
  {
    WALBERLA_ASSERT_NOT_NULLPTR( ptr );
  }

  template< typename DataType >
  DataType* get()
  {
    return static_cast< DataType* >( ptr_ );
  }

private:

  void *ptr_;

};

class ConstPrimitiveData : private NonCopyable
{
public:

  template< typename DataType >
  ConstPrimitiveData( const DataType* ptr ) : ptr_( ptr )
  {
    WALBERLA_ASSERT_NOT_NULLPTR( ptr );
  }

  template< typename DataType >
  const DataType* get()
  {
    return static_cast< const DataType* >( ptr_ );
  }

private:

  const void *ptr_;

};

}

class PrimitiveID;
class PrimitiveStorage;

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
  typedef internal::ConstPrimitiveData ConstPrimitiveData;

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
  void getNeighborVertices( std::vector< PrimitiveID > & neighborVertices ) const { neighborVertices.assign( neighborVertices_.begin(), neighborVertices_.end() ); }
  void getNeighborEdges   ( std::vector< PrimitiveID > & neighborEdges )    const { neighborEdges.assign   ( neighborEdges_.begin(),    neighborEdges_.end()    ); }
  void getNeighborFaces   ( std::vector< PrimitiveID > & neighborFaces )    const { neighborFaces.assign   ( neighborFaces_.begin(),    neighborFaces_.end()    ); }

  const std::vector< PrimitiveID > & neighborVertices() const { return neighborVertices_; }
  const std::vector< PrimitiveID > & neighborEdges()    const { return neighborEdges_; }
  const std::vector< PrimitiveID > & neighborFaces()    const { return neighborFaces_; }

  uint_t getNumNeighborVertices() const { return neighborVertices_.size(); }
  uint_t getNumNeighborEdges   () const { return neighborEdges_.size(); }
  uint_t getNumNeighborFaces   () const { return neighborFaces_.size(); }

  virtual void getLowerDimNeighbors ( std::vector< PrimitiveID > & lowerDimNeighbors )  const = 0;
  virtual void getHigherDimNeighbors( std::vector< PrimitiveID > & higherDimNeighbors ) const = 0;

  virtual const std::vector< PrimitiveID > & lowerDimNeighbors()  const = 0;
  virtual const std::vector< PrimitiveID > & higherDimNeighbors() const = 0;

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

  template< typename DataType, typename PrimitiveType >
  inline PrimitiveDataHandling< DataType, PrimitiveType >* getDataHandling( const PrimitiveDataID< DataType, PrimitiveType > & index ) const;

  /// Not public in order to guarantee that data is only added through the governing structure.
  /// This ensures valid DataIDs.
  template< typename DataType >
  inline void addData( const PrimitiveDataID< DataType, Primitive > & index,
	               const PrimitiveDataHandling< DataType, Primitive > & dataHandling );

  template< typename DataType, typename PrimitiveType >
  inline void genericAddData( const PrimitiveDataID< DataType, PrimitiveType > & index,
		              const PrimitiveDataHandling< DataType, PrimitiveType > & dataHandling,
		              const PrimitiveType * primitive );

  std::vector< PrimitiveID > neighborVertices_;
  std::vector< PrimitiveID > neighborEdges_;
  std::vector< PrimitiveID > neighborFaces_;

private:

  /// Holds a pointer to the actual data in the first entry and a pointer to the respective datahandling in the second entry.
  /// This way it is possible to loop over the data to for example serialize all registered data.
  std::vector< std::pair< PrimitiveData*, ConstPrimitiveData* > > data_;

  PrimitiveID primitiveID_;

};

// General methods for data and data handling retrieval
template< typename DataType, typename PrimitiveType >
DataType* Primitive::genericGetData( const PrimitiveDataID< DataType, PrimitiveType > & index ) const
{
  WALBERLA_ASSERT_LESS( index, data_.size(), "There is no data available for the specified index" );
  return data_[ index ].first->template get< DataType >();
}

template< typename DataType, typename PrimitiveType >
PrimitiveDataHandling< DataType, PrimitiveType >* Primitive::getDataHandling( const PrimitiveDataID< DataType, PrimitiveType > & index ) const
{
  WALBERLA_ASSERT_LESS( index, data_.size(), "There is no data handling available for the specified index" );
  return data_[ index ].second->template get< PrimitiveDataHandling< DataType, PrimitiveType > >();
}
///////////////////////////////////////////////////////

// Methods to retrieve data and data handling from primitives
template< typename DataType >
DataType* Primitive::getData( const PrimitiveDataID< DataType, Primitive > & index ) const
{
  return genericGetData< DataType >( index );
}

template< typename DataType >
PrimitiveDataHandling< DataType, Primitive >* Primitive::getDataHandling( const PrimitiveDataID< DataType, Primitive > & index ) const
{
  WALBERLA_ASSERT_LESS( index, data_.size(), "There is no data handling available for the specified index" );
  return data_[ index ].second->template get< PrimitiveDataHandling< DataType, Primitive > >();
}
/////////////////////////////////////////////////////////////

template< typename DataType >
void Primitive::addData( const PrimitiveDataID< DataType, Primitive > & index,
                         const PrimitiveDataHandling< DataType, Primitive > & dataHandling )
{
  genericAddData( index, dataHandling, this );
}

template< typename DataType, typename PrimitiveType >
void Primitive::genericAddData( const PrimitiveDataID< DataType, PrimitiveType > & index,
                                const PrimitiveDataHandling< DataType, PrimitiveType > & dataHandling,
	                        const PrimitiveType * const primitive )
{
  if( data_.size() <= index )
  {
    data_.resize( index + 1, std::make_pair< PrimitiveData*, ConstPrimitiveData* >( NULL, NULL ) );
  }

  WALBERLA_ASSERT_NOT_NULLPTR( primitive );

  data_[index].first = new PrimitiveData( dataHandling.initialize( primitive ) );
  data_[index].second = new ConstPrimitiveData( &dataHandling );
}


} // namespace hhg

