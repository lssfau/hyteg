
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "tinyhhg_core/tinyhhg.hpp"

#include "tinyhhg_core/indexing/VertexDoFIndexing.hpp"

namespace hhg {

/**
 * \page 04_Indexing Indexing functions and iterating over the domain
 *
 * \dontinclude tutorials/04_Indexing.cpp
 *
 * \brief In this tutorial we will iterate over simulation data
 *
 * \section intro Introduction
 *
 *
 * \section code Complete Program
 *
 * \include tutorials/04_Indexing.cpp
 *
 */
void IndexingTutorial()
{

  // The refinement level
  const uint_t level = 3;

  // The number of vertices on an edge on this level
  const uint_t numMicroverticesOnEdge = levelinfo::num_microvertices_per_edge( level );

  // While iterating over a face, the row size shrinks when increasing the row.
  // We need to keep track of it:
  uint_t innerRowSize = numMicroverticesOnEdge;

  // Iterating over the vertex DoFs (imagine a P1 FEM discretization) of a macro face
  for ( uint_t row = 0; row < numMicroverticesOnEdge; row++ )
  {
    for ( uint_t col = 0; col < innerRowSize; col++ )
    {
      // Calculate the array index using the apropriate index function.
      // For vertex DoFs (P1) on a macro face:
      const uint_t arrayIndex = indexing::vertexdof::macroface::index< level >( col, row );

      WALBERLA_LOG_INFO_ON_ROOT( "Array index for col = " << col << ", row = " << row << ": " << arrayIndex );
    }

    // Row size decreases in next row
    innerRowSize--;
  }

  // For stencil computations we typically only access the inner DoFs of a macro-primitive and
  // want to access their neighbors.

  // Convenient iteration using iterators:
  // Specify the sub-face to iterate on:
  const uint_t offsetToCenter = 1; // 0 for the whole face, 1 for the inner face (distance to border == 1), ...
  for ( const auto & it : indexing::vertexdof::macroface::Iterator< level, offsetToCenter >() )
  {
    // Get neighbors using the stencil versions of the indexing functions
    const uint_t vertexDoF      = indexing::vertexdof::macroface::indexFromVertex< level >( it.col(), it.row(), stencilDirection::VERTEX_C );
    const uint_t leftNeighbor   = indexing::vertexdof::macroface::indexFromVertex< level >( it.col(), it.row(), stencilDirection::VERTEX_W );
    const uint_t rightNeighbor  = indexing::vertexdof::macroface::indexFromVertex< level >( it.col(), it.row(), stencilDirection::VERTEX_E );

    WALBERLA_LOG_INFO_ON_ROOT( "Stencil array idx access: row = " << it.row() << ", col =  " << it.col() << ": " <<
                               "W = " << leftNeighbor << " C = " << vertexDoF << " E = " << rightNeighbor );
  }
}

}

int main( int argc, char** argv )
{
  walberla::mpi::Environment env( argc, argv );
  walberla::mpi::MPIManager::instance()->useWorldComm();
  hhg::IndexingTutorial();
  return 0;
}


