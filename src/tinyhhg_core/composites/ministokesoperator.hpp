#ifndef TINYHHG_MINISTOKESOPERATOR_HPP
#define TINYHHG_MINISTOKESOPERATOR_HPP

#include "tinyhhg_core/composites/ministokesfunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"

namespace hhg
{

class MiniStokesOperator
{
public:

  MiniStokesOperator(Mesh& mesh, size_t minLevel, size_t maxLevel)
    : A(mesh, minLevel, maxLevel),
      div_x(mesh, minLevel, maxLevel),
      div_y(mesh, minLevel, maxLevel),
      divT_x(mesh, minLevel, maxLevel),
      divT_y(mesh, minLevel, maxLevel)
  {
  }

  void apply(MiniStokesFunction& src, MiniStokesFunction& dst, size_t level, DoFType flag)
  {
    A.apply(src.u, dst.u, level, flag, Replace);
    divT_x.apply(src.p, dst.u, level, flag, Add);

    A.apply(src.v, dst.v, level, flag, Replace);
    divT_y.apply(src.p, dst.v, level, flag, Add);

    div_x.apply(src.u, dst.p, level, flag | DirichletBoundary, Replace);
    div_y.apply(src.v, dst.p, level, flag | DirichletBoundary, Add);
  }

  void save(const MiniStokesFunction& src, MiniStokesFunction& dst, std::ostream& out, size_t level, DoFType flag)
  {
    A.save(src.u, dst.u, out, level, flag);
    divT_x.save(src.p, dst.u, out, level, flag);

    A.save(src.v, dst.v, out, level, flag);
    divT_y.save(src.p, dst.v, out, level, flag);

    div_x.save(src.u, dst.p, out, level, flag);
    div_y.save(src.v, dst.p, out, level, flag);
  }

  P1BubbleLaplaceOperator A;
  P1BubbleToP1DivxOperator div_x;
  P1BubbleToP1DivyOperator div_y;
  P1ToP1BubbleDivTxOperator divT_x;
  P1ToP1BubbleDivTyOperator divT_y;
};

}

#endif //TINYHHG_MINISTOKESOPERATOR_HPP
