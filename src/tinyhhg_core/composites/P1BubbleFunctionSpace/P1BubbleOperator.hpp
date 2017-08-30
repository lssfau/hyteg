#pragma once

#include "tinyhhg_core/composites/P1BubbleFunctionSpace/P1BubbleFunction.hpp"

#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleOperator.hpp"

namespace hhg
{

template<class Opr_P1_to_P1, class Opr_B_to_B>
class P1BubbleOperator
{
public:

  P1BubbleOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : A_p1(storage, minLevel, maxLevel),
      A_b(storage, minLevel, maxLevel)
  {
  }

  void apply(P1BubbleFunction& src, P1BubbleFunction& dst, size_t level, DoFType flag, UpdateType updateType)
  {
//    WALBERLA_LOG_DEBUG_ON_ROOT("P1BubbleOperator::apply is missing off diagonal blocks!")

    A_p1.apply(src.p1, dst.p1, level, flag, updateType);
    A_b.apply(src.b, dst.b, level, flag, updateType);
  }

  void save(P1BubbleFunction& src, P1BubbleFunction& dst, Mat &mat, size_t level, DoFType flag)
  {
//    WALBERLA_LOG_DEBUG_ON_ROOT("P1BubbleOperator::apply is missing off diagonal blocks!")

    /*A_p1.save(src.p1, dst.p1, mat, level, flag); //TODO Implement
    A_b.save(src.b, dst.b, mat, level, flag);*/
  }

  Opr_P1_to_P1 A_p1;
  Opr_B_to_B   A_b;
};


typedef P1BubbleOperator<P1LaplaceOperator, BubbleLaplaceOperator> P1BubbleLaplaceOperator;

}
