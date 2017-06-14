
#pragma once

#include "core/debug/Debug.h"

namespace hhg {

using walberla::uint_t;
using walberla::blockforest::workload_t;
using walberla::blockforest::memory_t;

class SetupPrimitive : private walberla::NonCopyable
{
public:

  const PrimitiveID getPrimitiveID() const { return primitiveID_; }

        uint_t      getTargetRank()  const { return targetRank_; }
        void        setTargetRank( uint_t targetRank ) { targetRank_ = targetRank; }

        workload_t getWorkload() const { return workload_; }
        void       setWorkload( const workload_t w ) { WALBERLA_ASSERT_GREATER_EQUAL( w, static_cast< workload_t >(0) ); workload_ = w; }

        memory_t getMemory() const { return memory_; }
        void     setMemory( const memory_t m ) { WALBERLA_ASSERT_GREATER_EQUAL( m, static_cast< memory_t >(0) ); memory_ = m; }

protected:

  SetupPrimitive( const PrimitiveID & id ) :
    primitiveID_( id ), targetRank_( 0 ), workload_( 0 ), memory_( 0 )
  {}

private:

  PrimitiveID primitiveID_;
  uint_t      targetRank_;

  workload_t  workload_;
  memory_t    memory_;

};

}
