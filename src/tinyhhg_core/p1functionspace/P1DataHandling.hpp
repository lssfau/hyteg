
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "tinyhhg_core/p1functionspace/p1memory.hpp"

namespace hhg {

class VertexP1FunctionMemoryDataHandling : public OnlyInitializeDataHandling< VertexP1FunctionMemory, Vertex >
{
public:

  VertexP1FunctionMemory * initialize( const Vertex * const vertex ) const;

};

class EdgeP1FunctionMemoryDataHandling : public OnlyInitializeDataHandling< EdgeP1FunctionMemory, Edge >
{
public:

  EdgeP1FunctionMemory * initialize( const Edge * const edge ) const;

};

class FaceP1FunctionMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1FunctionMemory, Face >
{
public:

  FaceP1FunctionMemory * initialize( const Face * const face ) const;

};

}
