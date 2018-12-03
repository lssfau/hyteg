#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/ProlongationOperator.hpp"

namespace hhg {

class P1P1StokesToP1P1StokesProlongation : public ProlongationOperator< P1StokesFunction< real_t > >
{
public:

    typedef P1toP1LinearProlongation VelocityProlongation_T;
    typedef P1toP1LinearProlongation PressureProlongation_T;

    void prolongate ( const P1StokesFunction< real_t > & function, const uint_t & sourceLevel, const DoFType & flag ) override
    {
      prolongationOperator_( function.u, sourceLevel, flag );
      prolongationOperator_( function.v, sourceLevel, flag );
      prolongationOperator_( function.w, sourceLevel, flag );
      prolongationOperator_( function.p, sourceLevel, flag );
    }

private:

    P1toP1LinearProlongation prolongationOperator_;

};
}