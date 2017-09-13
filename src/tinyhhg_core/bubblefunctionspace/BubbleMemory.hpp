#pragma once

#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"

#include <string>


namespace hhg
{

template< typename ValueType >
class VertexBubbleFunctionMemory : public FunctionMemory< ValueType >
{
public:

  VertexBubbleFunctionMemory( const uint_t & numDependencies ) : FunctionMemory< ValueType >( numDependencies ) {}

  inline size_t getSize(size_t level) const
  {
    WALBERLA_UNUSED( level );
    return this->numDependencies_;
  }

};

template< typename ValueType >
class EdgeBubbleFunctionMemory : public FunctionMemory< ValueType >
{
public:

  EdgeBubbleFunctionMemory( const uint_t & numDependencies ) : FunctionMemory< ValueType >( numDependencies ) {}

  inline size_t getSize(size_t level) const
  {
    size_t num_cell_dofs = this->numDependencies_ * (2 * levelinfo::num_microedges_per_edge(level) - 1);
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
      data[level] = StencilStack{{std::unique_ptr< real_t[] >( new real_t[ getGrayStencilSize(level) ] ),
                                  std::unique_ptr< real_t[] >( new real_t[ getBlueStencilSize(level) ] )}};
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


template< typename ValueType >
class FaceBubbleFunctionMemory : public FunctionMemory< ValueType >
{
public:

  FaceBubbleFunctionMemory( const uint_t & numDependencies ) : FunctionMemory< ValueType >( numDependencies ) {}

  inline size_t getSize(size_t level) const
  {
    return levelinfo::num_microfaces_per_face(level);
  }

};

}
