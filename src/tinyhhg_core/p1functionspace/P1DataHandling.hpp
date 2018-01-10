
#pragma once

#include "tinyhhg_core/primitivedata/PrimitiveDataHandling.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"

namespace hhg {

template< typename ValueType >
class VertexP1StencilMemoryDataHandling : public StencilMemoryDataHandling< StencilMemory< ValueType >, Vertex >
{
 public:

  VertexP1StencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  inline std::shared_ptr< StencilMemory< ValueType > > initialize( const Vertex * const vertex ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class EdgeP1StencilMemoryDataHandling : public StencilMemoryDataHandling< StencilMemory< ValueType >, Edge >
{
 public:

  EdgeP1StencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  inline std::shared_ptr< StencilMemory< ValueType > > initialize( const Edge * const edge ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class FaceP1StencilMemoryDataHandling : public StencilMemoryDataHandling< StencilMemory< ValueType >, Face >
{
 public:

  FaceP1StencilMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  inline std::shared_ptr< StencilMemory< ValueType > > initialize( const Face * const face ) const;

 private:

  uint_t minLevel_;
  uint_t maxLevel_;

};



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


////////////////////
// Implementation //
////////////////////


template< typename ValueType >
std::shared_ptr< StencilMemory< ValueType > > VertexP1StencilMemoryDataHandling< ValueType >::initialize( const Vertex * const vertex ) const
{
  return std::make_shared< StencilMemory< ValueType > >( P1VertexStencilMemorySize, vertex->getNumNeighborEdges(), minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< StencilMemory< ValueType > > EdgeP1StencilMemoryDataHandling< ValueType >::initialize( const Edge * const edge ) const
{
  return std::make_shared< StencilMemory< ValueType > >( P1EdgeStencilMemorySize, edge->getNumNeighborFaces(), minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< StencilMemory< ValueType > > FaceP1StencilMemoryDataHandling< ValueType >::initialize( const Face * const ) const
{
  return std::make_shared< StencilMemory< ValueType > >( P1FaceStencilMemorySize, 0, minLevel_, maxLevel_ );
}


}
