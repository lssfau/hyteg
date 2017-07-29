#pragma once

#include "tinyhhg_core/support.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include <string>


namespace hhg
{

class VertexBubbleFunctionMemory
{
public:
  std::map<size_t, std::unique_ptr<real_t[]>> data;
  size_t num_deps_;

  inline std::unique_ptr<real_t[]>& addlevel(size_t level, size_t num_deps)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      this->num_deps_ = num_deps;
      data[level] = hhg::make_unique<real_t[]>(getSize(level));
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return num_deps_;
  }

};

class EdgeBubbleFunctionMemory
{
public:

  std::map<size_t, std::unique_ptr<real_t[]>> data;
  size_t num_deps_;

  inline std::unique_ptr<real_t[]>& addlevel(size_t level, size_t num_deps)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      this->num_deps_ = num_deps;
      data[level] = hhg::make_unique<real_t[]>(getSize(level));
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    size_t num_cell_dofs = num_deps_ * (2 * levelinfo::num_microedges_per_edge(level) - 1);
    return num_cell_dofs;
  }
};


class FaceBubbleStencilMemory
{
public:

  typedef std::array<std::unique_ptr<real_t[]>, 2> StencilStack;

  std::map<size_t, StencilStack> data;

  inline StencilStack& addlevel(size_t level)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      data[level] = StencilStack{{hhg::make_unique<real_t[]>(getGrayStencilSize(level)),
                                  hhg::make_unique<real_t[]>(getBlueStencilSize(level))}};
    }
    return data[level];
  }

  inline size_t getGrayStencilSize(size_t)
  {
    return 1;
  }

  inline size_t getBlueStencilSize(size_t)
  {
    return 1;
  }

};


class FaceBubbleFunctionMemory
{
public:

  std::map<size_t, std::unique_ptr<real_t[]>> data;

  inline std::unique_ptr<real_t[]>& addlevel(size_t level)
  {
    if (data.count(level)>0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      data[level] = hhg::make_unique<real_t[]>(getSize(level));
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return levelinfo::num_microfaces_per_face(level);
  }

};

}
