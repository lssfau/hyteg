#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "BubbleMemory.hpp"

namespace hhg {

template< typename ValueType >
class VertexBubbleFunctionMemoryDataHandling : public FunctionMemoryDataHandling< VertexBubbleFunctionMemory< ValueType >, Vertex >
{
public:

  VertexBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< VertexBubbleFunctionMemory< ValueType > > initialize( const Vertex * const vertex ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class EdgeBubbleFunctionMemoryDataHandling : public FunctionMemoryDataHandling< EdgeBubbleFunctionMemory< ValueType >, Edge >
{
public:

  EdgeBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< EdgeBubbleFunctionMemory< ValueType > > initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class FaceBubbleFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FaceBubbleFunctionMemory< ValueType >, Face >
{
public:

  FaceBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< FaceBubbleFunctionMemory< ValueType > > initialize( const Face * const face ) const;

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

template< typename ValueType >
std::shared_ptr< VertexBubbleFunctionMemory< ValueType > > VertexBubbleFunctionMemoryDataHandling< ValueType >::initialize( const Vertex * const vertex ) const
{
  return std::make_shared< VertexBubbleFunctionMemory< ValueType > >( bubbleVertexFunctionMemorySize, *vertex, minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< EdgeBubbleFunctionMemory< ValueType > > EdgeBubbleFunctionMemoryDataHandling< ValueType >::initialize( const Edge * const edge ) const
{
  return std::make_shared< EdgeBubbleFunctionMemory< ValueType > >( bubbleEdgeFunctionMemorySize, *edge, minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FaceBubbleFunctionMemory< ValueType > > FaceBubbleFunctionMemoryDataHandling< ValueType >::initialize( const Face * const face ) const
{
  return std::make_shared< FaceBubbleFunctionMemory< ValueType > >( bubbleFaceFunctionMemorySize, *face, minLevel_, maxLevel_ );
}

}
