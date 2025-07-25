/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Nils Kohl.
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

#include <typeinfo>

#include "hyteg/Levelinfo.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
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
template < typename VectorFunctionType >
class CSFVectorFunction;

template < typename FunctionTag_T, typename PrimitiveType >
inline uint_t numberOfInnerDoFs( const uint_t& level )
{
   // =============
   //  P1Functions
   // =============
   if constexpr ( std::is_same_v< FunctionTag_T, P1FunctionTag > )
   {
      // Vertex
      if constexpr ( std::is_same_v< PrimitiveType, Vertex > )
      {
         return 1;
      }
      // Edge
      else if constexpr ( std::is_same_v< PrimitiveType, Edge > )
      {
         return levelinfo::num_microvertices_per_edge( level ) - 2 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
      }
      // Face
      else if constexpr ( std::is_same_v< PrimitiveType, Face > )
      {
         return levelinfo::num_microvertices_per_face( level ) - 3 * numberOfInnerDoFs< P1FunctionTag, Edge >( level ) -
                3 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
      }
      // Cell
      else if constexpr ( std::is_same_v< PrimitiveType, Cell > )
      {
         return levelinfo::num_microvertices_per_cell( level ) - 4 * numberOfInnerDoFs< P1FunctionTag, Face >( level ) -
                6 * numberOfInnerDoFs< P1FunctionTag, Edge >( level ) - 4 * numberOfInnerDoFs< P1FunctionTag, Vertex >( level );
      }
   }

   // ==================
   //  EdgeDoFFunctions
   // ==================
   else if constexpr ( std::is_same_v< FunctionTag_T, EdgeDoFFunctionTag > )
   {
      // Vertex
      if constexpr ( std::is_same_v< PrimitiveType, Vertex > )
      {
         return 0;
      }
      // Edge
      else if constexpr ( std::is_same_v< PrimitiveType, Edge > )
      {
         return levelinfo::num_microedges_per_edge( level ) - 2 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
      }
      // Face
      else if constexpr ( std::is_same_v< PrimitiveType, Face > )
      {
         return levelinfo::num_microedges_per_face( level ) - 3 * numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) -
                3 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
      }
      // Cell
      else if constexpr ( std::is_same_v< PrimitiveType, Cell > )
      {
         return levelinfo::num_microedges_per_cell( level ) - 4 * numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level ) -
                6 * numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) -
                4 * numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
      }
   }

   // =============
   //  P2Functions
   // =============
   else if constexpr ( std::is_same_v< FunctionTag_T, P2FunctionTag > )
   {
      // Vertex
      if constexpr ( std::is_same_v< PrimitiveType, Vertex > )
      {
         return numberOfInnerDoFs< P1FunctionTag, Vertex >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level );
      }
      // Edge
      else if constexpr ( std::is_same_v< PrimitiveType, Edge > )
      {
         return numberOfInnerDoFs< P1FunctionTag, Edge >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level );
      }
      // Face
      else if constexpr ( std::is_same_v< PrimitiveType, Face > )
      {
         return numberOfInnerDoFs< P1FunctionTag, Face >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level );
      }
      // Cell
      else if constexpr ( std::is_same_v< PrimitiveType, Cell > )
      {
         return numberOfInnerDoFs< P1FunctionTag, Cell >( level ) + numberOfInnerDoFs< EdgeDoFFunctionTag, Cell >( level );
      }
   }

   // =============
   //  P0Functions
   // =============
   else if constexpr ( std::is_same_v< FunctionTag_T, P0FunctionTag > )
   {
      // Vertex
      if constexpr ( std::is_same_v< PrimitiveType, Vertex > )
      {
         return 0;
      }
      // Edge
      else if constexpr ( std::is_same_v< PrimitiveType, Edge > )
      {
         return 0;
      }
      // Face
      else if constexpr ( std::is_same_v< PrimitiveType, Face > )
      {
         return levelinfo::num_microfaces_per_face( level );
      }
      // Cell
      else if constexpr ( std::is_same_v< PrimitiveType, Cell > )
      {
         return levelinfo::num_microcells_per_cell( level );
      }
   }

   // ==============
   //  DG1Functions
   // ==============
   else if constexpr ( std::is_same_v< FunctionTag_T, DG1FunctionTag > )
   {
      // Vertex
      if constexpr ( std::is_same_v< PrimitiveType, Vertex > )
      {
         return 0;
      }
      // Edge
      else if constexpr ( std::is_same_v< PrimitiveType, Edge > )
      {
         return 0;
      }
      // Face
      else if constexpr ( std::is_same_v< PrimitiveType, Face > )
      {
         return 3 * levelinfo::num_microfaces_per_face( level );
      }
      // Cell
      else if constexpr ( std::is_same_v< PrimitiveType, Cell > )
      {
         return 4 * levelinfo::num_microcells_per_cell( level );
      }
   }

   // =======================
   //  P2PlusBubbleFunctions
   // =======================
   else if constexpr ( std::is_same_v< FunctionTag_T, P2PlusBubbleFunctionTag > )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "numberOfInnerDoFs() only gives a correct answer in 2D !!!" );
      return numberOfInnerDoFs< P2FunctionTag, PrimitiveType >( level ) +
             numberOfInnerDoFs< P0FunctionTag, PrimitiveType >( level );
   }

   // =========
   // NO MATCH
   // =========
   else
   {
      WALBERLA_ABORT( "Missing implementation in numberOfInnerDoFs(): typeID of Tag was '" << typeid( FunctionTag_T ).name()
                                                                                           << "'" );
   }

   return 0;
};

template < typename FunctionTag_T >
inline uint_t numberOfLocalDoFs( const PrimitiveStorage& primitiveStorage, const uint_t& level )
{
   // ============
   //  P0Function
   // ============
   if constexpr ( std::is_same_v< P0FunctionTag, FunctionTag_T > )
   {
      return primitiveStorage.hasGlobalCells() ?
                 numberOfInnerDoFs< P0FunctionTag, Cell >( level ) * primitiveStorage.getNumberOfLocalCells() :
                 numberOfInnerDoFs< P0FunctionTag, Face >( level ) * primitiveStorage.getNumberOfLocalFaces();
   }

   // =============
   //  DG1Function
   // =============
   if constexpr ( std::is_same_v< DG1FunctionTag, FunctionTag_T > )
   {
      return primitiveStorage.hasGlobalCells() ?
                 numberOfInnerDoFs< DG1FunctionTag, Cell >( level ) * primitiveStorage.getNumberOfLocalCells() :
                 numberOfInnerDoFs< DG1FunctionTag, Face >( level ) * primitiveStorage.getNumberOfLocalFaces();
   }

   // ============
   //  P1Function
   // ============
   else if constexpr ( std::is_same_v< P1FunctionTag, FunctionTag_T > )
   {
      return numberOfInnerDoFs< P1FunctionTag, Vertex >( level ) * primitiveStorage.getNumberOfLocalVertices() +
             numberOfInnerDoFs< P1FunctionTag, Edge >( level ) * primitiveStorage.getNumberOfLocalEdges() +
             numberOfInnerDoFs< P1FunctionTag, Face >( level ) * primitiveStorage.getNumberOfLocalFaces() +
             numberOfInnerDoFs< P1FunctionTag, Cell >( level ) * primitiveStorage.getNumberOfLocalCells();
   }

   // ==================
   //  P1VectorFunction
   // ==================
   else if constexpr ( std::is_same_v< P1VectorFunctionTag, FunctionTag_T > )
   {
      return ( primitiveStorage.hasGlobalCells() ? 3 : 2 ) * numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
   }

   // =================
   //  EdgeDoFFunction
   // =================
   else if constexpr ( std::is_same_v< EdgeDoFFunctionTag, FunctionTag_T > )
   {
      return numberOfInnerDoFs< EdgeDoFFunctionTag, Vertex >( level ) * primitiveStorage.getNumberOfLocalVertices() +
             numberOfInnerDoFs< EdgeDoFFunctionTag, Edge >( level ) * primitiveStorage.getNumberOfLocalEdges() +
             numberOfInnerDoFs< EdgeDoFFunctionTag, Face >( level ) * primitiveStorage.getNumberOfLocalFaces() +
             numberOfInnerDoFs< EdgeDoFFunctionTag, Cell >( level ) * primitiveStorage.getNumberOfLocalCells();
   }

   // ============
   //  P2Function
   // ============
   else if constexpr ( std::is_same_v< P2FunctionTag, FunctionTag_T > )
   {
      return numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level ) +
             numberOfLocalDoFs< EdgeDoFFunctionTag >( primitiveStorage, level );
   }

   // =======================
   //  P2PlusBubbleFunctions
   // =======================
   else if constexpr ( std::is_same_v< P2PlusBubbleFunctionTag, FunctionTag_T > )
   {
      if ( primitiveStorage.hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2PlusBubbleFunction restricted to 2D for the moment!" );
      }
      return numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level ) +
             numberOfLocalDoFs< P0FunctionTag >( primitiveStorage, level );
   }

   // ==================
   //  P2VectorFunction
   // ==================
   else if constexpr ( std::is_same_v< P2VectorFunctionTag, FunctionTag_T > )
   {
      return ( primitiveStorage.hasGlobalCells() ? 3 : 2 ) * numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level );
   }

   // =============================
   //  P2PlusBubbleVectorFunctions
   // =============================
   else if constexpr ( std::is_same_v< P2PlusBubbleVectorFunctionTag, FunctionTag_T > )
   {
      return ( primitiveStorage.hasGlobalCells() ? 3 : 2 ) *
             numberOfLocalDoFs< P2PlusBubbleFunctionTag >( primitiveStorage, level );
   }

   // ====================
   //  N1E1VectorFunction
   // ====================
   else if constexpr ( std::is_same_v< N1E1VectorFunctionTag, FunctionTag_T > )
   {
      return numberOfLocalDoFs< EdgeDoFFunctionTag >( primitiveStorage, level );
   }

   // ==================
   //  P1StokesFunction
   // ==================
   else if constexpr ( std::is_same_v< P1StokesFunctionTag, FunctionTag_T > )
   {
      return ( primitiveStorage.hasGlobalCells() ? 4 : 3 ) * numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
   }

   // ========================
   //  P2P1TaylorHoodFunction
   // ========================
   else if constexpr ( std::is_same_v< P2P1TaylorHoodFunctionTag, FunctionTag_T > )
   {
      return numberOfLocalDoFs< P2VectorFunctionTag >( primitiveStorage, level ) +
             numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
   }

   // =============================
   //  P2P1TaylorHoodBlockFunction
   // =============================
   else if constexpr ( std::is_same_v< P2P1TaylorHoodBlockFunctionTag, FunctionTag_T > )
   {
      return numberOfLocalDoFs< P2VectorFunctionTag >( primitiveStorage, level ) +
             numberOfLocalDoFs< P1FunctionTag >( primitiveStorage, level );
   }

   // ====================
   //  P2P2StokesFunction
   // ====================
   else if constexpr ( std::is_same_v< P2P2StokesFunctionTag, FunctionTag_T > )
   {
      return numberOfLocalDoFs< P2VectorFunctionTag >( primitiveStorage, level ) +
             numberOfLocalDoFs< P2FunctionTag >( primitiveStorage, level );
   }

   // ===================
   //  CCRStokesFunction
   // ===================
   else if constexpr ( std::is_same_v< CCRStokesFunctionTag, FunctionTag_T > )
   {
      return numberOfLocalDoFs< P2PlusBubbleVectorFunctionTag >( primitiveStorage, level ) +
             numberOfLocalDoFs< DG1FunctionTag >( primitiveStorage, level );
   }

   // =========
   // NO MATCH
   // =========
   else
   {
      WALBERLA_ABORT( "Missing implementation in numberOfLocalDoFs(): typeID of Tag was '" << typeid( FunctionTag_T ).name()
                                                                                           << "'" );
   }

   return 0;
};

/// variant of numberOfLocalDoFs that works on a function object
template < typename func_t >
inline uint_t numberOfLocalDoFs( const func_t& func, const uint_t& level )
{
   if constexpr ( std::is_base_of_v< BlockFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< GenericFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< dg::DGFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< CSFVectorFunction< dg::DGVectorFunction< typename func_t::valueType > >, func_t > ||
                  std::is_same_v< dg::DGVectorFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< P0Function< typename func_t::valueType >, func_t > ||
                  std::is_same_v< DG1Function< typename func_t::valueType >, func_t > ||
                  std::is_same_v< EGFunction< typename func_t::valueType >, func_t > )
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
                  std::is_same_v< CSFVectorFunction< dg::DGVectorFunction< typename func_t::valueType > >, func_t > ||
                  std::is_same_v< dg::DGVectorFunction< typename func_t::valueType >, func_t > ||
                  std::is_same_v< P0Function< typename func_t::valueType >, func_t > ||
                  std::is_same_v< DG1Function< typename func_t::valueType >, func_t > ||
                  std::is_same_v< EGFunction< typename func_t::valueType >, func_t > )
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

   // Logic for continous Galerkin type functions, for DG type functions
   // we do not need to subtract boundary dofs
   if constexpr ( !std::is_same_v< FunctionTag_T, P0FunctionTag > )
   {
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

/// Returns the polynomial degree of the basis functions of the corresponding function type, if implemented and meaningful.
template < typename FunctionType >
inline uint_t polynomialDegreeOfBasisFunctions()
{
   if ( std::is_same_v< typename FunctionType::Tag, P1FunctionTag > ||
        std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > )
   {
      return 1;
   }
   else if ( std::is_same_v< typename FunctionType::Tag, P2FunctionTag > ||
             std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag > )
   {
      return 2;
   }
   else
   {
      WALBERLA_ABORT( "polynomialDegreeOfBasisFunctions() not implemented or not meaningful for this function type." )
   }
}

/// \brief Prints the memory consumption of all functions with a specified scalar type.
///
/// Note 1: this function may be expensive since it performs global reduction (multiple times)
/// Note 2: must be called collectively by all threads
/// Note 3: might partly include memory allocated for stencils
///
template < typename ValueType = real_t >
void printFunctionAllocationInfo()
{
   // calculate actual allocated memory
   const double globalActualAllocatedMemory = double( FunctionMemory< ValueType >::getGlobalAllocatedMemoryInBytes() ) / 1e+09;
   const double minActualAllocatedMemory    = double( FunctionMemory< ValueType >::getMinLocalAllocatedMemoryInBytes() ) / 1e+09;
   const double maxActualAllocatedMemory    = double( FunctionMemory< ValueType >::getMaxLocalAllocatedMemoryInBytes() ) / 1e+09;

   WALBERLA_LOG_INFO_ON_ROOT( "====================== Function Allocation Info ======================" );
   WALBERLA_LOG_INFO_ON_ROOT( "The table entries are with respect to MPI processes." );
   WALBERLA_LOG_INFO_ON_ROOT( "Thus, min/max gives the minimal/maximal value of memory allocated on " );
   WALBERLA_LOG_INFO_ON_ROOT( "a process, while sum is the global amount over all processes." );
   WALBERLA_LOG_INFO_ON_ROOT( "=====================================================================" );
   WALBERLA_LOG_INFO_ON_ROOT( "                       +--------------+--------------+--------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "                       |          sum |          min |          max |" );
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------+--------------+--------------+--------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( " allocated memory (GB) | "
                              << std::setw( 12 ) << std::scientific << std::setprecision( 3 ) << globalActualAllocatedMemory
                              << " | " << std::setw( 12 ) << std::scientific << std::setprecision( 3 ) << minActualAllocatedMemory
                              << " | " << std::setw( 12 ) << std::scientific << std::setprecision( 3 ) << maxActualAllocatedMemory
                              << " | " );
   WALBERLA_LOG_INFO_ON_ROOT( " ----------------------+--------------+--------------+--------------+" );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

} // namespace hyteg
