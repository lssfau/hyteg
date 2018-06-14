
#pragma once

#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"

namespace hhg {

template< typename FunctionTag_T >
uint_t numberOfLocalDoFs( const PrimitiveStorage & primitiveStorage, const uint_t & level );

template<>
uint_t numberOfLocalDoFs< P1FunctionTag >( const PrimitiveStorage & primitiveStorage, const uint_t & level )
{
  const uint_t numDoFsOnSingleVertex = 1;
  const uint_t numDoFsOnSingleEdge   = levelinfo::num_microvertices_per_edge( level ) - 2 * numDoFsOnSingleVertex;
  const uint_t numDoFsOnSingleFace   = levelinfo::num_microvertices_per_face( level ) - 3 * numDoFsOnSingleEdge - 3 * numDoFsOnSingleVertex;
  const uint_t numDoFsOnSingleCell   = levelinfo::num_microvertices_per_cell( level ) - 4 * numDoFsOnSingleFace - 6 * numDoFsOnSingleEdge - 4 * numDoFsOnSingleVertex;
  return   numDoFsOnSingleVertex * primitiveStorage.getNumberOfLocalVertices()
           + numDoFsOnSingleEdge   * primitiveStorage.getNumberOfLocalEdges()
           + numDoFsOnSingleFace   * primitiveStorage.getNumberOfLocalFaces()
           + numDoFsOnSingleCell   * primitiveStorage.getNumberOfLocalCells();
}

template<>
uint_t numberOfLocalDoFs< EdgeDoFFunctionTag >( const PrimitiveStorage & primitiveStorage, const uint_t & level )
{
  const uint_t numDoFsOnSingleVertex = 0;
  const uint_t numDoFsOnSingleEdge   = levelinfo::num_microedges_per_edge( level ) - 2 * numDoFsOnSingleVertex;
  const uint_t numDoFsOnSingleFace   = levelinfo::num_microedges_per_face( level ) - 3 * numDoFsOnSingleEdge - 3 * numDoFsOnSingleVertex;
  return   numDoFsOnSingleVertex * primitiveStorage.getNumberOfLocalVertices()
           + numDoFsOnSingleEdge   * primitiveStorage.getNumberOfLocalEdges()
           + numDoFsOnSingleFace   * primitiveStorage.getNumberOfLocalFaces();
}

template<>
uint_t numberOfLocalDoFs< P2FunctionTag >( const PrimitiveStorage & primitiveStorage, const uint_t & level )
{
  return numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level ) + numberOfLocalDoFs< EdgeDoFFunctionTag >( primitiveStorage, level );
}

template<>
uint_t numberOfLocalDoFs< P1StokesFunctionTag >( const PrimitiveStorage & primitiveStorage, const uint_t & level )
{
  return 3 * numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
}

template<>
uint_t numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( const PrimitiveStorage & primitiveStorage, const uint_t & level )
{
  return 2 * numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level ) + numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
}

template< typename FunctionTag_T >
uint_t numberOfGlobalDoFs( const PrimitiveStorage & primitiveStorage, const uint_t & level )
{
  return walberla::mpi::allReduce( numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ), walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
}

}