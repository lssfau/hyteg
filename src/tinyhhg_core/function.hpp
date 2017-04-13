#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include "mesh.hpp"
#include "comm.hpp"
#include "operator.hpp"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"

#include <string>
#include <functional>

namespace hhg
{

class Function
{
public: 
  Function(const std::string& _name, Mesh& _mesh, size_t _minLevel, size_t _maxLevel)
    : name(_name), mesh(_mesh), minLevel(_minLevel), maxLevel(_maxLevel), memory_id(-1), rank(Comm::get().rk)
  {
  }

  virtual ~Function()
  {
  }

  virtual void interpolate(std::function<double(const hhg::Point3D&)>& expr, size_t level, size_t flag) = 0;

  template<class T, size_t N>
  void assign(std::array<double, N> scalars, std::array<T*, N> functions, size_t level, size_t flag);

  template<class T, size_t N>
  void add(std::array<double, N> scalars, std::array<T*, N> functions, size_t level, size_t flag);

  virtual double dot(Function& rhs, size_t level, size_t flag) = 0;

  virtual void apply(Operator& opr, Function& dst, size_t level, size_t flag) = 0;

  virtual void smooth_gs(Operator& opr, Function& rhs, size_t level, size_t flag) = 0;

  virtual void prolongate(size_t level, size_t flag) = 0;

  virtual void restrict(size_t level, size_t flag) = 0;

  std::string name;
  Mesh& mesh;
  size_t minLevel;
  size_t maxLevel;
  size_t memory_id;
  size_t rank;
};

}

#endif /* FUNCTION_HPP */