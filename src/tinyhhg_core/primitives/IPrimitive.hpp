
#pragma once

#include "domain_decomposition/IBlock.h"

namespace walberla {
namespace hhg {

class IPrimitive : private NonCopyable
{
public:

  typedef domain_decomposition::internal::BlockData PrimitiveData;

  inline const PrimitiveID & getId() { return primitiveID_; }

  bool operator==( const IBlock& rhs ) const;
  bool operator!=( const IBlock& rhs ) const { return !operator==( rhs ); }

  //**********************************************************************************************************************
  /*!
  *   Returns true if block data is allocated on this block at index "index"
  */
  //**********************************************************************************************************************
  bool isBlockDataAllocated( const ConstBlockDataID & index ) const;

  //**********************************************************************************************************************
  /*!
  *   For documentation of this function see "const T* IBlock::getData( const BlockDataID index ) const".
  */
  //**********************************************************************************************************************
  template< typename T >
  const T* getData( const ConstBlockDataID & index ) const;

  //**********************************************************************************************************************
  /*!
  *   Function for retrieving data - data that was previously added by the governing block storage data structure class
  *   (see class 'BlockStorage') via block data initialization functions.
  *   Data cannot be added manually by for example calling 'addData', but must be added through the governing block
  *   storage data structure that keeps track of assigning BlockDataIDs (BlockDataIDs are required for
  *   accessing/retrieving data stored within blocks).
  *   Retrieving data will fail if T is neither the same type nor a base class of the data that was added.
  */
  //**********************************************************************************************************************
  template< typename T >
  const T* getData( const BlockDataID & index ) const;

  //**********************************************************************************************************************
  /*!
  *   For documentation of this function see "const T* IBlock::getData( const BlockDataID index ) const".
  */
  //**********************************************************************************************************************
  template< typename T >
        T* getData( const BlockDataID & index );

  //**********************************************************************************************************************
  /*!
  *   Function for removing all data that corresponds to block data ID 'index'.
  *   Further calls to "getData" with 'index' will return NULL.
  */
  //**********************************************************************************************************************
  void deleteData( const BlockDataID & index );

protected:

  IPrimitive( PrimitiveStorage & storage, const PrimitiveID::IDType & id );

  virtual ~IPrimitive(); ///< Must not be made public! No one should be allowed to delete a variable of type 'IBlock*'

private:

  std::vector< PrimitiveData* > data_; ///< the data assigned to this block

  PrimitiveStorage & storage_; ///< governing block storage data structure

  PrimitiveID primitiveID_;

};



}
}
