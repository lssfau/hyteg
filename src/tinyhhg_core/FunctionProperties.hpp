
#pragma once

#include "tinyhhg_core/FunctionTraits.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg {

template < typename FunctionTag_T, typename PrimitiveType >
inline uint_t numberOfInnerDoFs( const uint_t& level );

template<>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Vertex >( const uint_t& )
{
   return 1;
}

template<>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Edge >( const uint_t& level )
{
   return levelinfo::num_microvertices_per_edge( level ) - 2 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
}

template<>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Face >( const uint_t& level )
{
   return levelinfo::num_microvertices_per_face( level ) - 3 * numberOfInnerDoFs< P1FunctionTag, Edge >( level ) - 3 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
}

template<>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Cell >( const uint_t& level )
{
   return levelinfo::num_microvertices_per_cell( level ) - 4 * numberOfInnerDoFs< P1FunctionTag, Face >( level ) -
   6 * numberOfInnerDoFs< P1FunctionTag, Edge >( level ) - 4 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
}

template<>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( const uint_t& )
{
   return 0;
}

template<>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( const uint_t & level )
{
   return levelinfo::num_microedges_per_edge( level ) - 2 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template<>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( const uint_t & level )
{
   return levelinfo::num_microedges_per_face( level ) - 3 * numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) - 3 *  numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template<>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Cell >( const uint_t & level )
{
   return levelinfo::num_microedges_per_cell( level ) - 4 * numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level )
          - 6 * numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) - 4 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template<>
inline uint_t numberOfInnerDoFs< P2FunctionTag, Vertex >( const uint_t& level )
{
   return numberOfInnerDoFs< P1FunctionTag, Vertex >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template<>
inline uint_t numberOfInnerDoFs< P2FunctionTag, Edge >( const uint_t& level )
{
  return numberOfInnerDoFs< P1FunctionTag, Edge >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level );
}

template<>
inline uint_t numberOfInnerDoFs< P2FunctionTag, Face >( const uint_t& level )
{
  return numberOfInnerDoFs< P1FunctionTag, Face >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level );
}

template<>
inline uint_t numberOfInnerDoFs< P2FunctionTag, Cell >( const uint_t& level )
{
  return numberOfInnerDoFs< P1FunctionTag, Cell >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Cell >( level );
}


template < typename FunctionTag_T >
inline uint_t numberOfLocalDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level );

template <>
inline uint_t numberOfLocalDoFs< P1FunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{

   return numberOfInnerDoFs< P1FunctionTag, Vertex >( level ) * primitiveStorage.getNumberOfLocalVertices() +
          numberOfInnerDoFs< P1FunctionTag, Edge >( level ) * primitiveStorage.getNumberOfLocalEdges() +
          numberOfInnerDoFs< P1FunctionTag, Face >( level ) * primitiveStorage.getNumberOfLocalFaces() +
          numberOfInnerDoFs< P1FunctionTag, Cell >( level ) * primitiveStorage.getNumberOfLocalCells();
}

template <>
inline uint_t numberOfLocalDoFs< EdgeDoFFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level ) * primitiveStorage.getNumberOfLocalVertices() +
          numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) * primitiveStorage.getNumberOfLocalEdges() +
          numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level ) * primitiveStorage.getNumberOfLocalFaces() +
          numberOfInnerDoFs< EdgeDoFFunctionTag, Cell >( level ) * primitiveStorage.getNumberOfLocalCells();
}

template <>
inline uint_t numberOfLocalDoFs< DGFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return levelinfo::num_microfaces_per_face( level ) * primitiveStorage.getNumberOfLocalFaces();
}

template <>
inline uint_t numberOfLocalDoFs< P2FunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level ) +
          numberOfLocalDoFs< EdgeDoFFunctionTag >( primitiveStorage, level );
}

template <>
inline uint_t numberOfLocalDoFs< P1StokesFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return ( primitiveStorage.hasGlobalCells() ? 4 : 3 ) * numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
}

template <>
inline uint_t numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return ( primitiveStorage.hasGlobalCells() ? 3 : 2 ) * numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level ) +
          numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
}

template <>
inline uint_t numberOfLocalDoFs< P2P2StokesFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return ( primitiveStorage.hasGlobalCells() ? 4 : 3 ) * numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level );
}

template < typename FunctionTag_T >
inline uint_t numberOfGlobalDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return walberla::mpi::allReduce( numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ),
                                    walberla::mpi::SUM,
                                    walberla::mpi::MPIManager::instance()->comm() );
}

} // namespace hhg