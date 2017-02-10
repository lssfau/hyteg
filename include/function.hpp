#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <string>
#include <functional>

#include "types/pointnd.hpp"
#include "types/flags.hpp"

#include "p1functionspace/p1operator.hpp"

namespace hhg
{

template<typename FunctionSpace>
class Function
{
public: 
  Function(const std::string& _name, FunctionSpace& _functionSpace, size_t _minLevel, size_t _maxLevel)
    : name(_name), functionSpace(_functionSpace), minLevel(_minLevel), maxLevel(_maxLevel)
  {
    memory_id = functionSpace.allocate(minLevel, maxLevel);
  }

  ~Function()
  {
    functionSpace.free(memory_id, minLevel, maxLevel);
  }

  void interpolate(std::function<double(const hhg::Point3D&)>& expr, size_t level, size_t flag = All)
  {
    functionSpace.interpolate(memory_id, expr, level, flag);
  }

  template<size_t N>
  void assign(std::array<double, N> scalars, std::array<Function*, N> functions, size_t level, size_t flag = All)
  {
    std::array<size_t, N> src_ids;
    for (size_t i = 0; i < N; ++i)
    {
      src_ids[i] = functions[i]->memory_id;
    }

    functionSpace.assign(scalars, src_ids, memory_id, level, flag);
  }

  template<size_t N>
  void add(std::array<double, N> scalars, std::array<Function*, N> functions, size_t level, size_t flag = All)
  {
    std::array<size_t, N> src_ids;
    for (size_t i = 0; i < N; ++i)
    {
      src_ids[i] = functions[i]->memory_id;
    }

    functionSpace.add(scalars, src_ids, memory_id, level, flag);
  }

  double dot(Function& rhs, size_t level, size_t flag = All)
  {
    return functionSpace.dot(memory_id, rhs.memory_id, level, flag);
  }

  void apply(P1LaplaceOperator& opr, Function& dst, size_t level, size_t flag = All)
  {
    functionSpace.apply(opr.id, memory_id, dst.memory_id, level, flag);
  }

  void smooth_gs(P1LaplaceOperator& opr, Function& rhs, size_t level, size_t flag = All)
  {
    functionSpace.smooth_gs(opr.id, memory_id, rhs.memory_id, level, flag);
  }

  void prolongate(size_t level, size_t flag = All)
  {
    functionSpace.prolongate(memory_id, level, flag);
  }

  void restrict(size_t level, size_t flag = All)
  {
    functionSpace.restrict(memory_id, level, flag);
  }

  std::string name;
  FunctionSpace& functionSpace;
  size_t minLevel;
  size_t maxLevel;
  size_t memory_id;
};

}

#endif /* FUNCTION_HPP */