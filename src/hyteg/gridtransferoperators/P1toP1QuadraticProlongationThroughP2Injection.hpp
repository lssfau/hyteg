
#pragma once

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/FunctionMemory.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"

namespace hyteg {

class P1toP1QuadraticProlongationThroughP2Injection : public ProlongationOperator< P1Function< real_t > >
{
public:
  P1toP1QuadraticProlongationThroughP2Injection( const std::shared_ptr< PrimitiveStorage >& storage,
                                                 const uint_t&                              p1MinLevel,
                                                 const uint_t&                              p1MaxLevel )
  : tmp_( "tmp", storage, p1MinLevel - 1, p1MaxLevel - 1)
  {
    WALBERLA_CHECK_GREATER( p1MinLevel, 2 );
  }

  inline void prolongate( const P1Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const override
  {
    tmp_.assign( function, sourceLevel - 1, All );
    quadraticProlongation_.prolongate( tmp_, sourceLevel - 1, flag );
    function.assign( tmp_, sourceLevel + 1, All );
  }

private:

    P2Function< real_t > tmp_;
    P2toP2QuadraticProlongation quadraticProlongation_;
};

}