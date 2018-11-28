
#pragma once

#include "tinyhhg_core/FunctionTraits.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

namespace hhg {

template < typename FunctionTag_T >
inline uint_t numberOfLocalDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level );

template <>
inline uint_t numberOfLocalDoFs< P1FunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   const uint_t numDoFsOnSingleVertex = 1;
   const uint_t numDoFsOnSingleEdge   = levelinfo::num_microvertices_per_edge( level ) - 2 * numDoFsOnSingleVertex;
   const uint_t numDoFsOnSingleFace =
       levelinfo::num_microvertices_per_face( level ) - 3 * numDoFsOnSingleEdge - 3 * numDoFsOnSingleVertex;
   const uint_t numDoFsOnSingleCell = levelinfo::num_microvertices_per_cell( level ) - 4 * numDoFsOnSingleFace -
                                      6 * numDoFsOnSingleEdge - 4 * numDoFsOnSingleVertex;
   return numDoFsOnSingleVertex * primitiveStorage.getNumberOfLocalVertices() +
          numDoFsOnSingleEdge * primitiveStorage.getNumberOfLocalEdges() +
          numDoFsOnSingleFace * primitiveStorage.getNumberOfLocalFaces() +
          numDoFsOnSingleCell * primitiveStorage.getNumberOfLocalCells();
}

template <>
inline uint_t numberOfLocalDoFs< EdgeDoFFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   const uint_t numDoFsOnSingleVertex = 0;
   const uint_t numDoFsOnSingleEdge   = levelinfo::num_microedges_per_edge( level ) - 2 * numDoFsOnSingleVertex;
   const uint_t numDoFsOnSingleFace =
       levelinfo::num_microedges_per_face( level ) - 3 * numDoFsOnSingleEdge - 3 * numDoFsOnSingleVertex;
   const uint_t numDoFsOnSingleCell =
       levelinfo::num_microedges_per_cell( level ) - 4 * numDoFsOnSingleFace
       - 6 * numDoFsOnSingleEdge - 4 * numDoFsOnSingleVertex;
   return numDoFsOnSingleVertex * primitiveStorage.getNumberOfLocalVertices() +
          numDoFsOnSingleEdge * primitiveStorage.getNumberOfLocalEdges() +
          numDoFsOnSingleFace * primitiveStorage.getNumberOfLocalFaces() +
          numDoFsOnSingleCell * primitiveStorage.getNumberOfLocalCells();
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
   return 3 * numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
}

template <>
inline uint_t numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return 2 * numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level ) +
          numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
}

template < typename FunctionTag_T >
inline uint_t numberOfGlobalDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return walberla::mpi::allReduce( numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ),
                                    walberla::mpi::SUM,
                                    walberla::mpi::MPIManager::instance()->comm() );
}

} // namespace hhg