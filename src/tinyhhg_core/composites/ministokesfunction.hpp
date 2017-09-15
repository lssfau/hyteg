#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/composites/P1BubbleFunctionSpace/P1BubbleFunction.hpp"

namespace hhg
{

template< typename ValueType >
class MiniStokesFunction
{
public:

  MiniStokesFunction(const std::string& _name, const std::shared_ptr< PrimitiveStorage > & storage, size_t _minLevel, size_t _maxLevel)
    : u(_name+"_u", storage, _minLevel, _maxLevel),
      v(_name+"_v", storage, _minLevel, _maxLevel),
      p(_name+"_p", storage, _minLevel, _maxLevel)
  {
  }

  void interpolate(std::function<real_t(const hhg::Point3D&)>& expr, size_t level, DoFType flag = All)
  {
    u.interpolate(expr, level, flag);
    v.interpolate(expr, level, flag);
    p.interpolate(expr, level, flag | DirichletBoundary);
  }

  void assign(const std::vector<walberla::real_t> scalars, const std::vector<MiniStokesFunction*> functions, size_t level, DoFType flag = All)
  {
    std::vector<P1BubbleFunction<real_t>*> functions_u;
    std::vector<P1BubbleFunction<real_t>*> functions_v;
    std::vector< P1Function< real_t > *> functions_p;

    for (auto& function : functions)
    {
      functions_u.push_back(&function->u);
      functions_v.push_back(&function->v);
      functions_p.push_back(&function->p);
    }

    u.assign(scalars, functions_u, level, flag);
    v.assign(scalars, functions_v, level, flag);
    p.assign(scalars, functions_p, level, flag | DirichletBoundary);
  }

  void add(const std::vector<walberla::real_t> scalars, const std::vector<MiniStokesFunction*> functions, size_t level, DoFType flag = All)
  {
    std::vector<P1BubbleFunction<real_t>*> functions_u;
    std::vector<P1BubbleFunction<real_t>*> functions_v;
    std::vector< P1Function< real_t > *> functions_p;

    for (auto& function : functions)
    {
      functions_u.push_back(&function->u);
      functions_v.push_back(&function->v);
      functions_p.push_back(&function->p);
    }

    u.add(scalars, functions_u, level, flag);
    v.add(scalars, functions_v, level, flag);
    p.add(scalars, functions_p, level, flag | DirichletBoundary);
  }

  walberla::real_t dot(MiniStokesFunction& rhs, size_t level, DoFType flag = All)
  {
    walberla::real_t sum = u.dot(rhs.u, level, flag);
    sum += v.dot(rhs.v, level, flag);
    sum += p.dot(rhs.p, level, flag | DirichletBoundary);
    return sum;
  }

  uint_t enumerate(size_t level, uint_t& num)
  {
    uint_t localSize = 0;
    localSize += u.enumerate(level, num);
    localSize += v.enumerate(level, num);
    localSize += p.enumerate(level, num);
    return localSize;
  }

//  void prolongate(size_t level, DoFType flag = All)
//  {
//    u.prolongate(level, flag);
//    v.prolongate(level, flag);
//    p.prolongate(level, flag | DirichletBoundary);
//  }
//
//  void restrict(size_t level, DoFType flag = All)
//  {
//    u.restrict(level, flag);
//    v.restrict(level, flag);
//    p.restrict(level, flag | DirichletBoundary);
//  }

  P1BubbleFunction< ValueType > u;
  P1BubbleFunction< ValueType > v;
  P1Function< ValueType > p;
};

}
