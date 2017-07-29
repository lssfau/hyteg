#pragma once

#include "tinyhhg_core/support.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include <string>


namespace hhg
{

//class VertexBubbleStencilMemory
//{
//public:
//
//  std::map<size_t, std::unique_ptr<real_t[]>> data;
//  size_t num_deps_;
//
//  inline std::unique_ptr<real_t[]>& addlevel(size_t level, size_t num_deps)
//  {
//    if (data.count(level)>0)
//      WALBERLA_LOG_WARNING("Level already exists.")
//    else
//    {
//      this->num_deps_ = num_deps;
//      data[level] = hhg::make_unique<real_t[]>(getSize(level));
//    }
//    return data[level];
//  }
//
//  inline size_t getSize(size_t level)
//  {
//    return levelinfo::num_microvertices_per_vertex(level) + num_deps_;
//  }
//
//};


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


//class EdgeBubbleStencilMemory
//{
//public:
//
//  std::map<size_t, std::unique_ptr<real_t[]>> data;
//
//  inline std::unique_ptr<real_t[]>& addlevel(size_t level)
//  {
//    //WALBERLA_LOG_DEVEL("EdgeP1StencilMemory, kind = " + std::to_string(this->type));
//    if (data.count(level)>0)
//      WALBERLA_LOG_WARNING("Level already exists.")
//    else
//    {
//      data[level] = hhg::make_unique<real_t[]>(getSize(level));
//    }
//    return data[level];
//  }
//
//  inline size_t getSize(size_t)
//  {
//    return 7;
//  }
//
//};


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

  inline size_t getSize(size_t)
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
