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

  std::shared_ptr< FunctionMemory< ValueType > > initialize( const Vertex * const vertex ) const override;

private:

uint_t minLevel_;
uint_t maxLevel_;

};


template< typename ValueType >
class EdgeDGFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory < ValueType >, Edge >
{
public:

  EdgeDGFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
    : minLevel_( minLevel ),
      maxLevel_( maxLevel )
  {}

  std::shared_ptr< FunctionMemory< ValueType > > initialize( const Edge * const vertex ) const override;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class FaceDGFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory < ValueType >, Face >
{
public:

  FaceDGFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
    : minLevel_( minLevel ),
      maxLevel_( maxLevel )
  {}

  std::shared_ptr<FunctionMemory<ValueType>> initialize( const Face *const face) const override;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > VertexDGFunctionMemoryDataHandling< ValueType >::initialize( const Vertex * const vertex ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( DGVertexFunctionMemorySize, vertex->getNumNeighborEdges(), minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > EdgeDGFunctionMemoryDataHandling< ValueType >::initialize( const Edge * const edge ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( DGEdgeFunctionMemorySize, edge->getNumNeighborFaces(), minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > FaceDGFunctionMemoryDataHandling< ValueType >::initialize( const Face * const face ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( DGFaceFunctionMemorySize, 0, minLevel_, maxLevel_ );
}


}