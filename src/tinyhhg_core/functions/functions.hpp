#ifndef TINYHHG_FUNCTIONS_HPP
#define TINYHHG_FUNCTIONS_HPP

#include "../function.hpp"
#include "../p1functionspace/p1function.hpp"

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

  real_t dot(P1StokesFunction& rhs, size_t level, size_t flag)
  {
    real_t sum = u.dot(rhs.u, level, flag);
    sum += v.dot(rhs.v, level, flag);
    sum += p.dot(rhs.p, level, flag);
    return sum;
  }

  void prolongate(size_t level, size_t flag)
  {
    u.prolongate(level, flag);
    v.prolongate(level, flag);
    p.prolongate(level, flag);
  }

  void restrict(size_t level, size_t flag)
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
