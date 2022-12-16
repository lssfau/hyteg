/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#pragma once

#include "hyteg/Levelinfo.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using namespace walberla::mpistubs;

// some forward declarations
template < typename value_t >
class BlockFunction;
template < typename value_t >
class GenericFunction;

template < typename FunctionTag_T, typename PrimitiveType >
inline uint_t numberOfInnerDoFs( const uint_t& level );

template <>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Vertex >( const uint_t& )
{
   return 1;
}

template <>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Edge >( const uint_t& level )
{
   return levelinfo::num_microvertices_per_edge( level ) - 2 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
}

template <>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Face >( const uint_t& level )
{
   return levelinfo::num_microvertices_per_face( level ) - 3 * numberOfInnerDoFs< P1FunctionTag, Edge >( level ) -
          3 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
}

template <>
inline uint_t numberOfInnerDoFs< P1FunctionTag, Cell >( const uint_t& level )
{
   return levelinfo::num_microvertices_per_cell( level ) - 4 * numberOfInnerDoFs< P1FunctionTag, Face >( level ) -
          6 * numberOfInnerDoFs< P1FunctionTag, Edge >( level ) - 4 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
}

template <>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( const uint_t& )
{
   return 0;
}

template <>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( const uint_t& level )
{
   return levelinfo::num_microedges_per_edge( level ) - 2 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template <>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( const uint_t& level )
{
   return levelinfo::num_microedges_per_face( level ) - 3 * numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) -
          3 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template <>
inline uint_t numberOfInnerDoFs< EdgeDoFFunctionTag, Cell >( const uint_t& level )
{
   return levelinfo::num_microedges_per_cell( level ) - 4 * numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level ) -
          6 * numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) -
          4 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template <>
inline uint_t numberOfInnerDoFs< P2FunctionTag, Vertex >( const uint_t& level )
{
   return numberOfInnerDoFs< P1FunctionTag, Vertex >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
}

template <>
inline uint_t numberOfInnerDoFs< P2FunctionTag, Edge >( const uint_t& level )
{
   return numberOfInnerDoFs< P1FunctionTag, Edge >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level );
}

template <>
inline uint_t numberOfInnerDoFs< P2FunctionTag, Face >( const uint_t& level )
{
   return numberOfInnerDoFs< P1FunctionTag, Face >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level );
}

template <>
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
inline uint_t numberOfLocalDoFs< P1VectorFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return ( primitiveStorage.hasGlobalCells() ? 3 : 2 ) * numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
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
inline uint_t numberOfLocalDoFs< FaceDoFFunction_old_Tag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
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
inline uint_t numberOfLocalDoFs< P2VectorFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return ( primitiveStorage.hasGlobalCells() ? 3 : 2 ) * numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level );
}

template <>
inline uint_t numberOfLocalDoFs< N1E1VectorFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return numberOfLocalDoFs< EdgeDoFFunctionTag >( primitiveStorage, level );
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
inline uint_t numberOfLocalDoFs< P2P1TaylorHoodBlockFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return numberOfLocalDoFs< P2VectorFunctionTag >( primitiveStorage, level ) +
          numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
}

template <>
inline uint_t numberOfLocalDoFs< P2P2StokesFunctionTag >( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return ( primitiveStorage.hasGlobalCells() ? 4 : 3 ) * numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level );
}

/// variant of numberOfLocalDoFs that works on a function object
template < typename func_t >
inline uint_t numberOfLocalDoFs( const func_t& func, const uint_t& level )
{
   if constexpr ( std::is_base_of_v< BlockFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< GenericFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< dg::DGFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< P0Function< typename func_t::valueType >, func_t > )
   {
      return func.getNumberOfLocalDoFs( level );
   }
   else
   {
      return numberOfLocalDoFs< typename func_t::Tag >( *( func.getStorage() ), level );
   }
}

template < typename FunctionTag_T >
inline uint_t numberOfGlobalDoFs( const PrimitiveStorage& primitiveStorage,
                                  const uint_t&           level,
                                  const MPI_Comm&         communicator = walberla::mpi::MPIManager::instance()->comm(),
                                  const bool&             onRootOnly   = false )
{
   if ( onRootOnly )
   {
      return walberla::mpi::reduce(
          numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ), walberla::mpi::SUM, 0, communicator );
   }
   else
   {
      return walberla::mpi::allReduce(
          numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ), walberla::mpi::SUM, communicator );
   }
}

/// variant of numberOfGlobalDoFs that works on a function object
template < typename func_t >
inline uint_t numberOfGlobalDoFs( const func_t&   func,
                                  const uint_t&   level,
                                  const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                  const bool&     onRootOnly   = false )
{
   if constexpr ( std::is_base_of_v< BlockFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< GenericFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< dg::DGFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< P0Function< typename func_t::valueType >, func_t > )
   {
      return func.getNumberOfGlobalDoFs( level, communicator, onRootOnly );
   }
   else
   {
      return numberOfGlobalDoFs< typename func_t::Tag >( *( func.getStorage() ), level, communicator, onRootOnly );
   }
}

template < typename FunctionTag_T >
inline uint_t minNumberOfLocalDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return walberla::mpi::allReduce( numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ),
                                    walberla::mpi::MIN,
                                    walberla::mpi::MPIManager::instance()->comm() );
}

template < typename FunctionTag_T >
inline uint_t maxNumberOfLocalDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return walberla::mpi::allReduce( numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ),
                                    walberla::mpi::MAX,
                                    walberla::mpi::MPIManager::instance()->comm() );
}

template < typename FunctionTag_T >
inline uint_t numberOfLocalInnerDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   uint_t boundaryPoints = 0;

   for ( const auto& it : primitiveStorage.getFaceIDs() )
   {
      if ( primitiveStorage.onBoundary( it, true ) )
      {
         boundaryPoints += numberOfInnerDoFs< FunctionTag_T, Face >( level );
      }
   }
   for ( const auto& it : primitiveStorage.getEdgeIDs() )
   {
      if ( primitiveStorage.onBoundary( it, true ) )
      {
         boundaryPoints += numberOfInnerDoFs< FunctionTag_T, Edge >( level );
      }
   }
   for ( const auto& it : primitiveStorage.getVertexIDs() )
   {
      if ( primitiveStorage.onBoundary( it, true ) )
      {
         boundaryPoints += numberOfInnerDoFs< FunctionTag_T, Vertex >( level );
      }
   }

   return numberOfLocalDoFs< FunctionTag_T >( primitiveStorage, level ) - boundaryPoints;
}

/**
 * calculates the global number of DoFs that are not on the domain boundary
 * performs a mpi all Reduce
 * @tparam FunctionTag_T type of Function (e.g. P1Function)
 * @param primitiveStorage
 * @param level refinement level
 * @return global number of points that are not on the domain boundary
 */
template < typename FunctionTag_T >
inline uint_t numberOfGlobalInnerDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return walberla::mpi::allReduce( numberOfLocalInnerDoFs< FunctionTag_T >( primitiveStorage, level ),
                                    walberla::mpi::SUM,
                                    walberla::mpi::MPIManager::instance()->comm() );
}

/**
 * calculates the min number of process-local DoFs over all processes that are not on the domain boundary
 * performs a mpi all Reduce
 * @tparam FunctionTag_T type of Function (e.g. P1Function)
 * @param primitiveStorage
 * @param level refinement level
 * @return process-local global minimum number of points that are not on the domain boundary
 */
template < typename FunctionTag_T >
inline uint_t minNumberOfLocalInnerDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return walberla::mpi::allReduce( numberOfLocalInnerDoFs< FunctionTag_T >( primitiveStorage, level ),
                                    walberla::mpi::MIN,
                                    walberla::mpi::MPIManager::instance()->comm() );
}

/**
 * calculates the max number of process-local DoFs over all processes that are not on the domain boundary
 * performs a mpi all Reduce
 * @tparam FunctionTag_T type of Function (e.g. P1Function)
 * @param primitiveStorage
 * @param level refinement level
 * @return process-local global maximum number of points that are not on the domain boundary
 */
template < typename FunctionTag_T >
inline uint_t maxNumberOfLocalInnerDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   return walberla::mpi::allReduce( numberOfLocalInnerDoFs< FunctionTag_T >( primitiveStorage, level ),
                                    walberla::mpi::MAX,
                                    walberla::mpi::MPIManager::instance()->comm() );
}

/// \brief Prints infos about the memory consumption of all functions.
///
/// Note 1: this function may be expensive since it performs global reduction (multiple times)
/// Note 2: must be called collectively by all threads
/// Note 3: currently restricted to functions of type real_t
/// Note 4: might partly include memory allocated for stencils
///
/// \param storage Primitive Storage containing the functions
/// \param verbosityLevel
///     0: only print memory consumption
///     1: additionally number of functions and DoFs of each type
///     2: additionally list names of allocated functions
void printFunctionAllocationInfo( const PrimitiveStorage& storage, const uint_t& verbosityLevel = 1 );

} // namespace hyteg
