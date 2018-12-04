
#pragma once

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/ProlongationOperator.hpp"

namespace hhg {

class P2P1StokesToP2P1StokesProlongation : public ProlongationOperator< P2P1TaylorHoodFunction< real_t > >
{
public:

    typedef P2toP2QuadraticProlongation VelocityProlongation_T;
    typedef P1toP1LinearProlongation    PressureProlongation_T;

    void prolongate ( const P2P1TaylorHoodFunction< real_t > & function, const uint_t & sourceLevel, const DoFType & flag ) override
    {
      quadraticProlongationOperator_.prolongate( function.u, sourceLevel, flag );
      quadraticProlongationOperator_.prolongate( function.v, sourceLevel, flag );
      linearProlongationOperator_.prolongate( function.p, sourceLevel, flag );
    }

private:

    P2toP2QuadraticProlongation quadraticProlongationOperator_;
    P1toP1LinearProlongation    linearProlongationOperator_;

};
}