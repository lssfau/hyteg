#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "BubbleMemory.hpp"

namespace hhg {

class VertexBubbleFunctionMemoryDataHandling : public OnlyInitializeDataHandling< VertexBubbleFunctionMemory, Vertex >
{
public:

  VertexBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  VertexBubbleFunctionMemory * initialize( const Vertex * const vertex ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class EdgeBubbleFunctionMemoryDataHandling : public OnlyInitializeDataHandling< EdgeBubbleFunctionMemory, Edge >
{
public:

  EdgeBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  EdgeBubbleFunctionMemory * initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceBubbleFunctionMemoryDataHandling : public OnlyInitializeDataHandling< FaceBubbleFunctionMemory, Face >
{
public:

  FaceBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  FaceBubbleFunctionMemory * initialize( const Face * const face ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceBubbleStencilMemoryDataHandling : public OnlyInitializeDataHandling< FaceBubbleStencilMemory, Face >
{
 public:

  FaceBubbleStencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  FaceBubbleStencilMemory * initialize( const Face * const face ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

}
