
#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/RestrictionOperator.hpp"

namespace hhg {

class P1P1StokesToP1P1StokesRestriction : RestrictionOperator< P1StokesFunction< real_t > >
{
public:

    void restrict ( const P1StokesFunction< real_t > & function, const uint_t & sourceLevel, const DoFType & flag )
    {
      restrictionOperator_.restrict( function.u, sourceLevel, flag );
      restrictionOperator_.restrict( function.v, sourceLevel, flag );
      restrictionOperator_.restrict( function.w, sourceLevel, flag );
      restrictionOperator_.restrict( function.p, sourceLevel, flag );
    }

private:

    P1toP1LinearRestriction restrictionOperator_;

};
}