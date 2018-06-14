
#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"

namespace hhg {

class P1P1StokesToP1P1StokesProlongation
{
public:

    void operator() ( const P1StokesFunction< real_t > & function, const uint_t & sourceLevel, const DoFType & flag )
    {
      prolongationOperator_( function.u, sourceLevel, flag );
      prolongationOperator_( function.v, sourceLevel, flag );
      prolongationOperator_( function.p, sourceLevel, flag | DirichletBoundary );
    }

private:

    P1toP1LinearProlongation prolongationOperator_;

};
}