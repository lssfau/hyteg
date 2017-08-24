#pragma once

#include "tinyhhg_core/support.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include <string>


namespace hhg
{

class VertexP1StencilMemory
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
    return levelinfo::num_microvertices_per_vertex(level) + num_deps_;
  }

};


class VertexP1FunctionMemory : public FunctionMemory
{
public:

  VertexP1FunctionMemory( const uint_t & numDependencies ) : FunctionMemory( numDependencies ) {}

  inline size_t getSize(size_t level) const
  {
    return levelinfo::num_microvertices_per_vertex(level) + numDependencies_;
  }

};


class EdgeP1StencilMemory
{
public:

  std::map<size_t, std::unique_ptr<real_t[]>> data;

  inline std::unique_ptr<real_t[]>& addlevel(size_t level)
  {
    //WALBERLA_LOG_DEVEL("EdgeP1StencilMemory, kind = " + std::to_string(this->type));
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
    return 7;
  }

};


class EdgeP1FunctionMemory : public FunctionMemory
{
public:

  EdgeP1FunctionMemory( const uint_t & numDependencies ) : FunctionMemory( numDependencies ) {}

  inline size_t getSize(size_t level) const
  {
    size_t num_dofs_per_edge = levelinfo::num_microvertices_per_edge(level);
    return num_dofs_per_edge + numDependencies_*(num_dofs_per_edge-1);
  }
};


class FaceP1StencilMemory
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
    return 7;
  }

};


class FaceP1FunctionMemory : public FunctionMemory
{
public:

  FaceP1FunctionMemory( const uint_t & numDependencies ) : FunctionMemory( numDependencies ) {}

  inline size_t getSize(size_t level) const
  {
    return levelinfo::num_microvertices_per_face(level);
  }

};

}
