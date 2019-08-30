#pragma once

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/StokesOperatorTraits.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P1StokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
public:

  typedef P1ConstantLaplaceOperator VelocityOperator_T;
  typedef P1ConstantLaplaceOperator PressureOperator_T;

  P1StokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
  : Operator( storage, minLevel, maxLevel )
  , A( storage, minLevel, maxLevel )
  , div_x( storage, minLevel, maxLevel )
  , div_y( storage, minLevel, maxLevel )
  , div_z( storage, minLevel, maxLevel )
  , divT_x( storage, minLevel, maxLevel )
  , divT_y( storage, minLevel, maxLevel )
  , divT_z( storage, minLevel, maxLevel )
  , pspg( storage, minLevel, maxLevel )
  , hasGlobalCells_( storage->hasGlobalCells() )
  {}

  void apply(const P1StokesFunction<real_t>& src, const P1StokesFunction<real_t>& dst,const size_t level, DoFType flag) const
  {
    A.apply(src.u, dst.u, level, flag, Replace);
    divT_x.apply(src.p, dst.u, level, flag, Add);

    A.apply(src.v, dst.v, level, flag, Replace);
    divT_y.apply(src.p, dst.v, level, flag, Add);

    if ( hasGlobalCells_ )
    {
      A.apply(src.w, dst.w, level, flag, Replace);
      divT_z.apply(src.p, dst.w, level, flag, Add);
    }

    div_x.apply(src.u, dst.p, level, flag, Replace);
    div_y.apply(src.v, dst.p, level, flag, Add);

    if ( hasGlobalCells_ )
    {
      div_z.apply(src.w, dst.p, level, flag, Add);
    }

    pspg.apply(src.p, dst.p, level, flag, Add);
  }

  P1ConstantLaplaceOperator A;
  P1DivxOperator div_x;
  P1DivyOperator div_y;
  P1DivzOperator div_z;
  P1DivTxOperator divT_x;
  P1DivTyOperator divT_y;
  P1DivTzOperator divT_z;
  P1PSPGOperator pspg;
  bool hasGlobalCells_;
};

template<>
struct has_pspg_block< P1StokesOperator > {
    static const bool value = true;
};

}
