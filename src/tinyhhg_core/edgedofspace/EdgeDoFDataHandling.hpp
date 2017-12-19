#pragma once

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/primitives/all.hpp"
#include "tinyhhg_core/levelinfo.hpp"



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

}/// namespace hhg
