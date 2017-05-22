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

  template<size_t Level>
  void interpolate(std::function<double(const hhg::Point3D&)>& expr, DoFType flag = All)
  {
    u.interpolate<Level>(expr, flag);
    v.interpolate<Level>(expr, flag);
    p.interpolate<Level>(expr, flag | DirichletBoundary);
  }

  template<size_t Level>
  void assign(const std::vector<walberla::real_t> scalars, const std::vector<P1StokesFunction*> functions, DoFType flag = All)
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

    u.assign<Level>(scalars, functions_u, flag);
    v.assign<Level>(scalars, functions_v, flag);
    p.assign<Level>(scalars, functions_p, flag | DirichletBoundary);
  }

  template<size_t Level>
  void add(const std::vector<walberla::real_t> scalars, const std::vector<P1StokesFunction*> functions, DoFType flag = All)
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

    u.add<Level>(scalars, functions_u, flag);
    v.add<Level>(scalars, functions_v, flag);
    p.add<Level>(scalars, functions_p, flag | DirichletBoundary);
  }

  template<size_t Level>
  walberla::real_t dot(P1StokesFunction& rhs, DoFType flag = All)
  {
    walberla::real_t sum = u.dot<Level>(rhs.u, flag);
    sum += v.dot<Level>(rhs.v, flag);
    sum += p.dot<Level>(rhs.p, flag | DirichletBoundary);
    return sum;
  }

  template<size_t Level>
  void prolongate(DoFType flag = All)
  {
    u.prolongate<Level>(flag);
    v.prolongate<Level>(flag);
    p.prolongate<Level>(flag | DirichletBoundary);
  }

  template<size_t Level>
  void restrict(DoFType flag = All)
  {
    u.restrict<Level>(flag);
    v.restrict<Level>(flag);
    p.restrict<Level>(flag | DirichletBoundary);
  }

  P1Function u;
  P1Function v;
  P1Function p;
};

}

#endif //TINYHHG_FUNCTIONS_HPP
