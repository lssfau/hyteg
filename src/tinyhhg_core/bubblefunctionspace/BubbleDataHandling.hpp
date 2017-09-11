#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "BubbleMemory.hpp"

namespace hhg {

class VertexBubbleFunctionMemoryDataHandling : public FunctionMemoryDataHandling< VertexBubbleFunctionMemory< real_t >, Vertex >
{
public:

  VertexBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< VertexBubbleFunctionMemory< real_t > > initialize( const Vertex * const vertex ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class EdgeBubbleFunctionMemoryDataHandling : public FunctionMemoryDataHandling< EdgeBubbleFunctionMemory< real_t >, Edge >
{
public:

  EdgeBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< EdgeBubbleFunctionMemory< real_t > > initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceBubbleFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FaceBubbleFunctionMemory< real_t >, Face >
{
public:

  FaceBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< FaceBubbleFunctionMemory< real_t > > initialize( const Face * const face ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceBubbleStencilMemoryDataHandling : public OnlyInitializeDataHandling< FaceBubbleStencilMemory, Face >
{
 public:

  FaceBubbleStencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< FaceBubbleStencilMemory > initialize( const Face * const face ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

}
