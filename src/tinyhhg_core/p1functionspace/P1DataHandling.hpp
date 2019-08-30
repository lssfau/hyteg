
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

namespace hyteg {

class VertexP1LocalMatrixMemoryDataHandling : public OnlyInitializeDataHandling< VertexP1LocalMatrixMemory, Vertex >
{
public:

VertexP1LocalMatrixMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

std::shared_ptr< VertexP1LocalMatrixMemory > initialize( const Vertex * const vertex ) const;

private:

uint_t minLevel_;
uint_t maxLevel_;

};

class EdgeP1LocalMatrixMemoryDataHandling : public OnlyInitializeDataHandling< EdgeP1LocalMatrixMemory, Edge >
{
public:

  EdgeP1LocalMatrixMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< EdgeP1LocalMatrixMemory > initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceP1LocalMatrixMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1LocalMatrixMemory, Face >
{
public:

  FaceP1LocalMatrixMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< FaceP1LocalMatrixMemory > initialize( const Face * const face ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceP1PolynomialMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1PolynomialMemory, Face >
{
public:

  FaceP1PolynomialMemoryDataHandling( const uint_t & maxDegree ) : maxDegree_( maxDegree ) {}

  std::shared_ptr< FaceP1PolynomialMemory > initialize( const Face * const face ) const;

private:
  uint_t maxDegree_;

};


}
