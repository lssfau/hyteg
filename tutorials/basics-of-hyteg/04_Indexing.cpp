/*
 * Copyright (c) 2017-2019 Dominik Bartuschat, Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "hyteg/Levelinfo.hpp"

namespace hyteg {

/**
 * \page 04_Indexing Tutorial 04 - Indexing and iterating
 *
 * \dontinclude tutorials/basics-of-hyteg/04_Indexing.cpp
 *
 * \brief In this tutorial we will iterate over simulation data
 *
 * \section T04-Indexing-intro Introduction
 *
 * To allow for convenient iteration over the degrees of freedom of a macro-primitive in a FEM simulation
 * we introduce indexing functions that translate logical coordinates of the DoFs to array indices.
 *
 * These functions exist for the different types of DoFs, such as vertex (for P1 functions),
 * edge (used together with the vertex DoFs for P2) or cell (for DG functions in 3D), in combination with different
 * types of macro-primitives.
 *
 * In this tutorial we will show different approaches to iterate over the vertex DoFs of a macro face.
 *
 * \section T04-Indexing-forloop A Simple Loop
 *
 * First we need to specify the refinement level, which defines the number of DoFs in the macro primitive.
 *
 * \snippet tutorials/basics-of-hyteg/04_Indexing.cpp Level
 *
 * The logical layouts of the DoFs for the different primitive types are defined in the documentation.
 * For vertex DoFs on a macro face with refinement level 3, the DoFs are numbered as follows:
 *
 * \verbatim

   ( column, row )

   ( 0, 8 )
   ( 0, 7 ) ( 1, 7 )

     ...      ...
     ...      ...       ...

    ( 0, 1 ) ( 1, 1 ) ( 2, 1 ) ( 3, 1 ) ( 4, 1 ) ( 5, 1 ) ( 6, 1 ) ( 7, 1 )
    ( 0, 0 ) ( 1, 0 ) ( 2, 0 ) ( 3, 0 ) ( 4, 0 ) ( 5, 0 ) ( 6, 0 ) ( 7, 0 ) ( 8, 0 )

   \endverbatim
 *
 * Therefore, to build a loop over the triangle structure, we need to calculate the number of DoFs
 * in the bottom row i.e. the width (and height) of the structure:
 *
 * \snippet tutorials/basics-of-hyteg/04_Indexing.cpp LevelInfo
 *
 * The loop over the logical indices could be implemented as follows:
 *
 * \snippet tutorials/basics-of-hyteg/04_Indexing.cpp FirstLoop
 *
 * In the loop body, the logial coordinates are translated to the index in the corresponding array
 * using the index function for vertex DoFs on a macro-face:
 *
 * \snippet tutorials/basics-of-hyteg/04_Indexing.cpp IndexFunction
 *
 * \section T04-Indexing-iterators Iterators
 *
 * However, often we do not want to iterate over the complete macro-face but rather over the
 * inner DoFs or over the borders of the macro-face in order to collect data that shall be communicated
 * to neighbor primitives.
 *
 * Iterators can be used to avoid having to write complex loops by hand.
 *
 * For example to iterate over the inner macro-face, we can call the iterator
 * for vertex DoFs on a macro-face with the parameter offsetToCenter set to 1:
 *
 * \snippet tutorials/basics-of-hyteg/04_Indexing.cpp Iterator
 *
 * The Index variable \code it \endcode contains the current (logical) coordinates of the iterator in
 * the macro face.
 *
 * For convenient access to the neighboring DoFs of a certain DoF, there are indexing functions that return
 * the array index of neighboring DoFs of a certain direction:
 *
 * \snippet tutorials/basics-of-hyteg/04_Indexing.cpp Stencil
 *
 * \section T04-Indexing-code Complete Program
 *
 * \include tutorials/basics-of-hyteg/04_Indexing.cpp
 *
 */
void IndexingTutorial()
{

  /// [Level]
  const uint_t level = 3;
  /// [Level]

  /// [LevelInfo]
  const uint_t numMicroverticesOnEdge = levelinfo::num_microvertices_per_edge( level );
  /// [LevelInfo]

  /// [FirstLoop]
  uint_t innerRowSize = numMicroverticesOnEdge;

  for ( idx_t row = 0; row < idx_t( numMicroverticesOnEdge ); row++ )
  {
    for ( idx_t col = 0; col < idx_t( innerRowSize ); col++ )
    {
      // Calculate the array index using the appropriate index function.
      // For vertex DoFs (P1) on a macro face:
      /// [IndexFunction]
      const uint_t arrayIndex = vertexdof::macroface::index( level, col, row );
      /// [IndexFunction]

      WALBERLA_LOG_INFO_ON_ROOT( "Array index for col = " << col << ", row = " << row << ": " << arrayIndex );
    }

    // Row size decreases in next row
    innerRowSize--;
  }
  /// [FirstLoop]

  /// [Iterator]
  const uint_t offsetToCenter = 1; // 0 for the whole face, 1 for the inner face (distance to border == 1), ...

  for ( const auto & it : vertexdof::macroface::Iterator( level, offsetToCenter ) )
  /// [Iterator]
  {
    /// [Stencil]
    // Get neighbors using the stencil versions of the indexing functions
    const uint_t vertexDoF      = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(),
                                                                                  stencilDirection::VERTEX_C );
    const uint_t leftNeighbor   = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(),
                                                                                  stencilDirection::VERTEX_W );
    const uint_t rightNeighbor  = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(),
                                                                                  stencilDirection::VERTEX_E );
    /// [Stencil]

    WALBERLA_LOG_INFO_ON_ROOT( "Stencil array idx access: row = " << it.y() << ", col =  " << it.x() << ": " <<
                               "W = " << leftNeighbor << " C = " << vertexDoF << " E = " << rightNeighbor );
  }

}

}

int main( int argc, char** argv )
{
  walberla::mpi::Environment env( argc, argv );
  walberla::mpi::MPIManager::instance()->useWorldComm();
  hyteg::IndexingTutorial();
  return 0;
}


