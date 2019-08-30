
#pragma once

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"

namespace hyteg {

class P2P1StokesToP2P1StokesProlongation : public ProlongationOperator< P2P1TaylorHoodFunction< real_t > >
{
public:

    typedef P2toP2QuadraticProlongation VelocityProlongation_T;
    typedef P1toP1LinearProlongation    PressureProlongation_T;

    void prolongate ( const P2P1TaylorHoodFunction< real_t > & function, const uint_t & sourceLevel, const DoFType & flag ) const override
    {
      quadraticProlongationOperator_.prolongate( function.u, sourceLevel, flag );
      quadraticProlongationOperator_.prolongate( function.v, sourceLevel, flag );
      quadraticProlongationOperator_.prolongate( function.w, sourceLevel, flag );
      linearProlongationOperator_.prolongate( function.p, sourceLevel, flag );
    }

private:

    P2toP2QuadraticProlongation quadraticProlongationOperator_;
    P1toP1LinearProlongation    linearProlongationOperator_;

};
}