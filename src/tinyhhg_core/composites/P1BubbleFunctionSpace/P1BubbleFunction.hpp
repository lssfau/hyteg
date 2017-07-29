#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/bubblefunctionspace/BubbleFunction.hpp"

namespace hhg
{

class P1BubbleFunction
{
public:

  P1BubbleFunction(const std::string& _name, const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : p1(_name+"_p1", storage, minLevel, maxLevel),
      b(_name+"_b", storage, minLevel, maxLevel)
  {
  }

  void interpolate(std::function<real_t(const hhg::Point3D&)>& expr, size_t level, DoFType flag = All)
  {
    p1.interpolate(expr, level, flag);
//    b.interpolate(expr, level, flag);
  }

  void assign(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All)
  {
    std::vector<P1Function*> functions_p1;
    std::vector<BubbleFunction*> functions_b;

    for (auto& function : functions)
    {
      functions_p1.push_back(&function->p1);
      functions_b.push_back(&function->b);
    }

    p1.assign(scalars, functions_p1, level, flag);
    b.assign(scalars, functions_b, level, flag);
  }

  void add(const std::vector<walberla::real_t> scalars, const std::vector<P1BubbleFunction*> functions, size_t level, DoFType flag = All)
  {
    std::vector<P1Function*> functions_p1;
    std::vector<BubbleFunction*> functions_b;

    for (auto& function : functions)
    {
      functions_p1.push_back(&function->p1);
      functions_b.push_back(&function->b);
    }

    p1.add(scalars, functions_p1, level, flag);
    b.add(scalars, functions_b, level, flag);
  }

  walberla::real_t dot(P1BubbleFunction& rhs, size_t level, DoFType flag = All)
  {
    walberla::real_t sum = p1.dot(rhs.p1, level, flag);
    sum += b.dot(rhs.b, level, flag);
    return sum;
  }

  void prolongate(size_t level, DoFType flag = All)
  {
    WALBERLA_ASSERT(false, "P1BubbleFunction::prolongate is not implemented!");
  }

  void restrict(size_t level, DoFType flag = All)
  {
    WALBERLA_ASSERT(false, "P1BubbleFunction::restrict is not implemented!");
  }

  P1Function p1;
  BubbleFunction b;
};

}
