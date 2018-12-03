
#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/RestrictionOperator.hpp"

namespace hhg {

class P1P1StokesToP1P1StokesRestriction : RestrictionOperator< P1StokesFunction< real_t > >
{
public:

    typedef P1toP1LinearRestriction VelocityRestriction_T;
    typedef P1toP1LinearRestriction PressureRestriction_T;

    void restrict ( const P1StokesFunction< real_t > & function, const uint_t & sourceLevel, const DoFType & flag )
    {
      restrictionOperator_( function.u, sourceLevel, flag );
      restrictionOperator_( function.v, sourceLevel, flag );
      restrictionOperator_( function.w, sourceLevel, flag );
      restrictionOperator_( function.p, sourceLevel, flag );
    }

private:

    P1toP1LinearRestriction restrictionOperator_;

};
}