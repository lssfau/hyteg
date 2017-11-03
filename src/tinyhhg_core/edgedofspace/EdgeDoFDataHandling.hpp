#pragma once

#include "tinyhhg_core/edgedofspace/EdgeDoFMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

namespace hhg {

template< typename ValueType >
class EdgeDoFMacroVertexFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory< ValueType >, Vertex >
{
public:

  EdgeDoFMacroVertexFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  inline std::shared_ptr< FunctionMemory< ValueType > > initialize( const Vertex * const vertex ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class EdgeDoFMacroEdgeFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory< ValueType >, Edge >
{
public:

  EdgeDoFMacroEdgeFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  inline std::shared_ptr< FunctionMemory< ValueType > > initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class EdgeDoFMacroFaceFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory< ValueType >, Face >
{
public:

  EdgeDoFMacroFaceFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  inline std::shared_ptr< FunctionMemory< ValueType > > initialize( const Face * const face ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

////////////////////
// Implementation //
////////////////////

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > EdgeDoFMacroVertexFunctionMemoryDataHandling< ValueType >::initialize( const Vertex * const vertex ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( EdgeDoFMacroVertexFunctionMemorySize, vertex->getNumNeighborEdges(), minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > EdgeDoFMacroEdgeFunctionMemoryDataHandling< ValueType >::initialize( const Edge * const edge ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( EdgeDoFMacroEdgeFunctionMemorySize, edge->getNumNeighborFaces(), minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > EdgeDoFMacroFaceFunctionMemoryDataHandling< ValueType >::initialize( const Face * const ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( EdgeDoFMacroFaceFunctionMemorySize, 0, minLevel_, maxLevel_ );
}

}
