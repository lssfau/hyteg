
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg {

class VertexP1FunctionMemoryDataHandling : public OnlyInitializeDataHandling< VertexP1FunctionMemory, Vertex >
{
public:

  VertexP1FunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  VertexP1FunctionMemory * initialize( const Vertex * const vertex ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class EdgeP1FunctionMemoryDataHandling : public OnlyInitializeDataHandling< EdgeP1FunctionMemory, Edge >
{
public:

  EdgeP1FunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  EdgeP1FunctionMemory * initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceP1FunctionMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1FunctionMemory, Face >
{
public:

  FaceP1FunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  FaceP1FunctionMemory * initialize( const Face * const face ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

}
