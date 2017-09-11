
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

namespace hhg {

class VertexP1FunctionMemoryDataHandling : public FunctionMemoryDataHandling< VertexP1FunctionMemory< real_t >, Vertex >
{
public:

  VertexP1FunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< VertexP1FunctionMemory< real_t > > initialize( const Vertex * const vertex ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class EdgeP1FunctionMemoryDataHandling : public FunctionMemoryDataHandling< EdgeP1FunctionMemory< real_t >, Edge >
{
public:

  EdgeP1FunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< EdgeP1FunctionMemory< real_t > > initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceP1FunctionMemoryDataHandling : public FunctionMemoryDataHandling< FaceP1FunctionMemory< real_t >, Face >
{
public:

  FaceP1FunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< FaceP1FunctionMemory< real_t > > initialize( const Face * const face ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};


class VertexP1StencilMemoryDataHandling : public OnlyInitializeDataHandling< VertexP1StencilMemory, Vertex >
{
 public:

  VertexP1StencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< VertexP1StencilMemory > initialize( const Vertex * const vertex ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class EdgeP1StencilMemoryDataHandling : public OnlyInitializeDataHandling< EdgeP1StencilMemory, Edge >
{
 public:

  EdgeP1StencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< EdgeP1StencilMemory > initialize( const Edge * const edge ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceP1StencilMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1StencilMemory, Face >
{
 public:

  FaceP1StencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< FaceP1StencilMemory > initialize( const Face * const face ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

}
