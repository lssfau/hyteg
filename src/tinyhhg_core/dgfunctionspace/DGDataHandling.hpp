#pragma once

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"



namespace hhg{

template< typename ValueType >
class VertexDGFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory < ValueType >, Vertex >
{
public:

VertexBubbleFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
  : minLevel_( minLevel ),
    maxLevel_( maxLevel )
{}

std::shared_ptr< FunctionMemory< ValueType > > initialize( const Vertex * const vertex ) const;

private:

uint_t minLevel_;
uint_t maxLevel_;

};

template< typename ValueType >
std::shared_ptr< VertexBubbleFunctionMemory< ValueType > > VertexBubbleFunctionMemoryDataHandling< ValueType >::initialize( const Vertex * const vertex ) const
{
  return std::make_shared< VertexBubbleFunctionMemory< ValueType > >( bubbleVertexFunctionMemorySize, vertex->getNumNeighborFaces(), minLevel_, maxLevel_ );
}



}