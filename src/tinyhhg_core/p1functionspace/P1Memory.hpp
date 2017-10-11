#pragma once

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"

#include "tinyhhg_core/levelinfo.hpp"

#include <string>


namespace hhg
{

inline uint_t P1VertexFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  return levelinfo::num_microvertices_per_vertex( level ) + numDependencies;
}

inline uint_t P1EdgeFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  size_t num_dofs_per_edge = levelinfo::num_microvertices_per_edge( level );
  return num_dofs_per_edge + numDependencies * ( num_dofs_per_edge - 1 );
}

inline uint_t P1FaceFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  return levelinfo::num_microvertices_per_face(level);
}

template< typename ValueType >
using VertexP1FunctionMemory = FunctionMemory< ValueType >;

template< typename ValueType >
using EdgeP1FunctionMemory = FunctionMemory< ValueType >;

template< typename ValueType >
using FaceP1FunctionMemory = FunctionMemory< ValueType >;

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
      data[level] = std::unique_ptr< real_t[] >( new real_t[ getSize(level) ] );
    }
    return data[level];
  }

  inline size_t getSize(size_t level)
  {
    return levelinfo::num_microvertices_per_vertex(level) + num_deps_;
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
      data[level] = std::unique_ptr< real_t[] >( new real_t[ getSize(level) ] );
    }
    return data[level];
  }

  inline size_t getSize(size_t)
  {
    return 7;
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
      data[level] = std::unique_ptr< real_t[] >( new real_t[ getSize(level) ] );
    }
    return data[level];
  }

  inline size_t getSize(size_t)
  {
    return 7;
  }

};


}
