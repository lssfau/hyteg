
#pragma once

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"

namespace hhg {

class P2P1StokesToP2P1StokesRestriction
{
public:

    void operator() ( const P2P1TaylorHoodFunction< real_t > & function, const uint_t & sourceLevel, const DoFType & flag )
    {
      quadraticRestrictionOperator_( function.u, sourceLevel, flag );
      quadraticRestrictionOperator_( function.v, sourceLevel, flag );
      linearRestrictionOperator_   ( function.p, sourceLevel, flag | DirichletBoundary );
    }

private:

    P2toP2QuadraticRestriction quadraticRestrictionOperator_;
    P1toP1LinearRestriction    linearRestrictionOperator_;

};
}