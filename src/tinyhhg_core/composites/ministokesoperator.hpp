#ifndef TINYHHG_MINISTOKESOPERATOR_HPP
#define TINYHHG_MINISTOKESOPERATOR_HPP

#include "tinyhhg_core/composites/ministokesfunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"

namespace hhg
{

class MiniStokesOperator
{
public:

  MiniStokesOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : A(storage, minLevel, maxLevel),
      div_x_p1(storage, minLevel, maxLevel),
      div_y_p1(storage, minLevel, maxLevel),
      div_x_b(storage, minLevel, maxLevel),
      div_y_b(storage, minLevel, maxLevel),
      divT_x_p1(storage, minLevel, maxLevel),
      divT_y_p1(storage, minLevel, maxLevel),
      divT_x_b(storage, minLevel, maxLevel),
      divT_y_b(storage, minLevel, maxLevel)
  {
  }

  void apply(MiniStokesFunction& src, MiniStokesFunction& dst, size_t level, DoFType flag)
  {
    A.apply(src.u, dst.u, level, flag, Replace);
    divT_x_p1.apply(src.p, dst.u.p1, level, flag, Add);
    divT_x_b.apply(src.p, dst.u.b, level, flag, Add);

    A.apply(src.v, dst.v, level, flag, Replace);
    divT_y_p1.apply(src.p, dst.v.p1, level, flag, Add);
    divT_y_b.apply(src.p, dst.v.b, level, flag, Add);

    div_x_p1.apply(src.u.p1, dst.p, level, flag | DirichletBoundary, Replace);
    div_x_b.apply(src.u.b, dst.p, level, flag | DirichletBoundary, Add);
    div_y_p1.apply(src.v.p1, dst.p, level, flag | DirichletBoundary, Add);
    div_y_b.apply(src.v.b, dst.p, level, flag | DirichletBoundary, Add);
  }

//  void save(const MiniStokesFunction& src, MiniStokesFunction& dst, std::ostream& out, size_t level, DoFType flag)
//  {
//    A.save(src.u, dst.u, out, level, flag);
//    divT_x.save(src.p, dst.u, out, level, flag);
//
//    A.save(src.v, dst.v, out, level, flag);
//    divT_y.save(src.p, dst.v, out, level, flag);
//
//    div_x.save(src.u, dst.p, out, level, flag);
//    div_y.save(src.v, dst.p, out, level, flag);
//  }

  P1BubbleLaplaceOperator A;

  BubbleToP1DivTxOperator div_x_b;
  BubbleToP1DivTyOperator div_y_b;
  P1DivxOperator div_x_p1;
  P1DivyOperator div_y_p1;


  P1ToBubbleDivxOperator divT_x_b;
  P1ToBubbleDivxOperator divT_y_b;
  P1DivTxOperator divT_x_p1;
  P1DivTyOperator divT_y_p1;

};

}

#endif //TINYHHG_MINISTOKESOPERATOR_HPP
