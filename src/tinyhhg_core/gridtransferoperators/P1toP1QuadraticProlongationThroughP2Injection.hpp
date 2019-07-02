
#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/gridtransferoperators/ProlongationOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"

namespace hhg {

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