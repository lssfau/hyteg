
#pragma once

#include "core/debug/Debug.h"

namespace hhg {

using walberla::uint_t;
using walberla::blockforest::workload_t;
using walberla::blockforest::memory_t;

class SetupPrimitiveStorage;

class SetupPrimitive : private walberla::NonCopyable
{
public:

  virtual ~SetupPrimitive() {}

  const PrimitiveID getPrimitiveID() const { return primitiveID_; }

        workload_t getWorkload() const { return workload_; }
        void       setWorkload( const workload_t w ) { WALBERLA_ASSERT_GREATER_EQUAL( w, static_cast< workload_t >(0) ); workload_ = w; }

        memory_t getMemory() const { return memory_; }
        void     setMemory( const memory_t m ) { WALBERLA_ASSERT_GREATER_EQUAL( m, static_cast< memory_t >(0) ); memory_ = m; }

  /// @name Neighborhood
  /// Access to neighbors of either lower or higher dimension.
  ///@{
  PrimitiveID::const_iterator beginLowerDimNeighbors() const { return lowerDimNeighbors_.begin(); }
  PrimitiveID::const_iterator endLowerDimNeighbors()   const { return lowerDimNeighbors_.end(); }

  PrimitiveID::const_iterator beginHigherDimNeighbors() const { return higherDimNeighbors_.begin(); }
  PrimitiveID::const_iterator endHigherDimNeighbors()   const { return higherDimNeighbors_.end(); }

  uint_t getNumLowerDimNeighbors() const  { return lowerDimNeighbors_.size(); }
  uint_t getNumHigherDimNeighbors() const { return higherDimNeighbors_.size(); }
  ///@}

protected:

  SetupPrimitive( const PrimitiveID & id ) :
    primitiveID_( id ), workload_( 0 ), memory_( 0 )
  {}

  std::vector< PrimitiveID > lowerDimNeighbors_;
  std::vector< PrimitiveID > higherDimNeighbors_;

private:

  PrimitiveID primitiveID_;

  workload_t  workload_;
  memory_t    memory_;

};

}
