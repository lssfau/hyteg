#ifndef TINYHHG_FUNCTIONS_HPP
#define TINYHHG_FUNCTIONS_HPP

#include "tinyhhg_core/p1functionspace/p1function.hpp"

namespace hhg
{

class P1StokesFunction
{
public:

  P1StokesFunction(const std::string& _name, Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : u(_name+"_u", _mesh, _minLevel, _maxLevel),
      v(_name+"_v", _mesh, _minLevel, _maxLevel),
      p(_name+"_p", _mesh, _minLevel, _maxLevel)
  {
  }

  void assign(const std::vector<real_t> scalars, const std::vector<P1StokesFunction*> functions, size_t level, size_t flag = All)
  {
    std::vector<P1Function*> functions_u;
    std::vector<P1Function*> functions_v;
    std::vector<P1Function*> functions_p;

    for (auto& function : functions)
    {
      functions_u.push_back(&function->u);
      functions_v.push_back(&function->v);
      functions_p.push_back(&function->p);
    }

    u.assign(scalars, functions_u, level, flag);
    v.assign(scalars, functions_v, level, flag);
    p.assign(scalars, functions_p, level, flag);
  }

  void add(const std::vector<real_t> scalars, const std::vector<P1StokesFunction*> functions, size_t level, size_t flag = All)
  {
    std::vector<P1Function*> functions_u;
    std::vector<P1Function*> functions_v;
    std::vector<P1Function*> functions_p;

    for (auto& function : functions)
    {
      functions_u.push_back(&function->u);
      functions_v.push_back(&function->v);
      functions_p.push_back(&function->p);
    }

    u.add(scalars, functions_u, level, flag);
    v.add(scalars, functions_v, level, flag);
    p.add(scalars, functions_p, level, flag);
  }

  real_t dot(P1StokesFunction& rhs, size_t level, size_t flag = All)
  {
    real_t sum = u.dot(rhs.u, level, flag);
    sum += v.dot(rhs.v, level, flag);
    sum += p.dot(rhs.p, level, flag);
    return sum;
  }

  void prolongate(size_t level, size_t flag = All)
  {
    u.prolongate(level, flag);
    v.prolongate(level, flag);
    p.prolongate(level, flag);
  }

  void restrict(size_t level, size_t flag = All)
  {
    u.restrict(level, flag);
    v.restrict(level, flag);
    p.restrict(level, flag);
  }

  P1Function u;
  P1Function v;
  P1Function p;
};

}

#endif //TINYHHG_FUNCTIONS_HPP
