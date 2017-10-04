#pragma once

#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"

#include <string>


namespace hhg
{

inline uint_t bubbleVertexFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( level );
  return numDependencies;
}

inline uint_t bubbleEdgeFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  size_t num_cell_dofs = numDependencies * ( 2 * levelinfo::num_microedges_per_edge( level ) - 1 );
  return num_cell_dofs;
}

inline uint_t bubbleFaceFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( numDependencies );
  return levelinfo::num_microfaces_per_face( level );
}

template< typename ValueType >
using VertexBubbleFunctionMemory = FunctionMemory< ValueType >;

template< typename ValueType >
using EdgeBubbleFunctionMemory = FunctionMemory< ValueType >;

template< typename ValueType >
using FaceBubbleFunctionMemory = FunctionMemory< ValueType >;


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


}
