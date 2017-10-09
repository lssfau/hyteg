#pragma once

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "DGMemory.hpp"

using namespace hhg;

namespace hhg{

template< typename ValueType >
class VertexDGFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory < ValueType >, Vertex >
{
public:

VertexDGFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
  : minLevel_( minLevel ),
    maxLevel_( maxLevel )
{}

std::shared_ptr< FunctionMemory< ValueType > > initialize( const Vertex * const vertex ) const;

private:

uint_t minLevel_;
uint_t maxLevel_;

};

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > VertexDGFunctionMemoryDataHandling< ValueType >::initialize( const Vertex * const vertex ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( DGVertexFunctionMemorySize, vertex->getNumNeighborFaces(), minLevel_, maxLevel_ );
}



}