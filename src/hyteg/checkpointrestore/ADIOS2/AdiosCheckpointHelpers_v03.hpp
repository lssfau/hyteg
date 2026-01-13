/*
 * Copyright (c) 2023-2025 Marcus Mohr.
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

#include <adios2.h>

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/types/Concepts.hpp"

namespace hyteg::adiosCheckpointHelpers_v03 {

/// Generate an ADIOS2 variable name for a function
///
/// The name is composed by combining the name of the function with the refinement
/// level for which data are to be exported.
std::string generateVariableName( const std::string& functionName, uint_t level )
{
   std::stringstream sStream;
   sStream << functionName << "-LVL" << level;
   return sStream.str();
}

/// Computes the maximum data size of a certain function on either a face or cell primitive including halos
template < concepts::concrete_primitive primitive_t, template < typename > class func_t, concepts::value_type value_t >
   requires concepts::fe_function_scalar< func_t< value_t > > &&
            (std::is_same_v< primitive_t, Face > || std::is_same_v< primitive_t, Cell >)
uint_t getMaxDataBlockSize( uint_t level )
{
   uint_t blockSize = 0u;

   // Treat case of scalar P1Function
   if constexpr ( std::is_same_v< func_t< value_t >, vertexdof::VertexDoFFunction< value_t > > )
   {
      if constexpr ( std::is_same_v< primitive_t, Face > )
      {
         blockSize = levelinfo::num_microvertices_per_face( level );
      }
      else
      {
         blockSize = levelinfo::num_microvertices_per_cell( level );
      }
   }

   // Treat case of scalar EdgeDoFFunction
   else if constexpr ( std::is_same_v< func_t< value_t >, EdgeDoFFunction< value_t > > )
   {
      if constexpr ( std::is_same_v< primitive_t, Face > )
      {
         blockSize =
             3u * ( ( ( levelinfo::num_microedges_per_edge( level ) + 1u ) * levelinfo::num_microedges_per_edge( level ) ) / 2u );
      }
      else
      {
         blockSize = 6u * ( levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) ) +
                     ( levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) - 1u ) );
      }
   }

   // Function not supported, yet
   else
   {
      WALBERLA_ABORT( "Implement blocksize computation for another kind of fe function!" );
   }

   return blockSize;
}

template < template < typename > class func_t, concepts::value_type value_t, typename action_t >
void doSomethingForAFunctionOnAllHighestDimensionalPrimitives( adios2::IO&                            io,
                                                               adios2::Engine&                        engine,
                                                               const func_t< value_t >&               function,
                                                               uint_t                                 minLevel,
                                                               uint_t                                 maxLevel,
                                                               const std::map< PrimitiveID, uint_t >& rowIndices,
                                                               action_t                               action,
                                                               uint_t                                 step = 0U )
{
   const auto& storage = function.getStorage();

   if ( storage->hasGlobalCells() )
   {
      for ( const auto& macroIterator : storage->getCells() )
      {
         const Cell& cell = *macroIterator.second;
         for ( uint_t level = minLevel; level <= maxLevel; ++level )
         {
            action( io, engine, function, level, rowIndices.at( cell.getID() ), cell, step );
         }
      }
   }

   else
   {
      for ( const auto& macroIterator : storage->getFaces() )
      {
         const Face& face = *macroIterator.second;
         for ( uint_t level = minLevel; level <= maxLevel; ++level )
         {
            action( io, engine, function, level, rowIndices.at( face.getID() ), face, step );
         }
      }
   }
}

/// Define an ADIOS2 variable for a scalar FEFunction
///
/// We use a gloabl array for this. The number of rows will be the global number of highest dimensional
/// primitives in the storage on which the function lives, i.e. faces or cells. The number of columns
/// will be the number of DoFs allocated for the function on a single such primitive on the given
/// refinement level including the halo.
template < template < typename > class func_t, concepts::value_type value_t >
   requires concepts::fe_function_scalar< func_t< value_t > >
void generateVariableForScalarFunction( adios2::IO&              io,
                                        const func_t< value_t >& function,
                                        uint_t                   level,
                                        uint_t                   globalNumberOfHighestDimPrimitives )
{
   // create a name for the variable
   std::string varName = generateVariableName( function.getFunctionName(), level );

   // determine the size of data we have (depends on the function, the level and the type of primitive)
   // and also the number of highest-dimensional primitives we own locally
   uint_t      chunkSize     = 0u;
   uint_t      numPrimitives = 0u;
   const auto& storage       = function.getStorage();
   if ( storage->hasGlobalCells() )
   {
      numPrimitives = storage->getNumberOfLocalCells();
      chunkSize     = getMaxDataBlockSize< Cell, func_t, value_t >( level );
   }
   else
   {
      numPrimitives = storage->getNumberOfLocalFaces();
      chunkSize     = getMaxDataBlockSize< Face, func_t, value_t >( level );
   }

   // if we do not have any local primitives, we will not put anything into the global array
   uint_t localSize = numPrimitives > 0 ? chunkSize : 0u;

   WALBERLA_LOG_PROGRESS( "Going to define variable '"
                          << varName << "' with\nshape .... {" << globalNumberOfHighestDimPrimitives << ", " << chunkSize << "}\n"
                          << "start .... {0, 0}\ncount .... {" << numPrimitives << ", " << localSize << "}" );

   // now we can define the variable; note that we will provide the local offsets later with SetSelection()
   io.DefineVariable< value_t >(
       varName, { globalNumberOfHighestDimPrimitives, chunkSize }, { 0, 0 }, { numPrimitives, localSize } );
};

/// Schedule data associated with a primitive and a scalar FE function for export
template < template < typename > class func_t, concepts::value_type value_t >
   requires concepts::fe_function_scalar< func_t< value_t > >
void exportVariableForScalarFunction( adios2::IO&              io,
                                      adios2::Engine&          engine,
                                      const func_t< value_t >& function,
                                      uint_t                   level,
                                      uint_t                   rowIndexInAdiosArray,
                                      const Primitive&         primitive,
                                      uint_t                   step = 0U )
{
   WALBERLA_UNUSED( step );

   // need a name for the global array variable
   std::string varName = generateVariableName( function.getFunctionName(), level );

   // check that associated variable exists in IO object
   adios2::Variable< value_t > varDoFData = io.InquireVariable< value_t >( varName );
   if ( !varDoFData )
   {
      WALBERLA_ABORT( "ADIOS2 variable '" << varName << "' does not exist!" );
   }

   switch ( primitive.getType() )
   {
   case Primitive::FACE: {
      const Face& face = dynamic_cast< const Face& >( primitive );

      // figure out what the size of our data chunk actually is
      uint_t size    = face.getData( function.getFaceDataID() )->getSize( level );
      uint_t maxSize = getMaxDataBlockSize< Face, func_t, value_t >( level );
      WALBERLA_CHECK_LESS_EQUAL( size, maxSize );

      // now we can set the position for this primitive's data chunk inside the global array
      varDoFData.SetSelection( { { rowIndexInAdiosArray, 0 }, { 1, size } } );

      // schedule the data for exporting
      engine.Put( varDoFData, face.getData( function.getFaceDataID() )->getPointer( level ) );

      break;
   }
   case Primitive::CELL: {
      const Cell& cell = dynamic_cast< const Cell& >( primitive );

      // figure out what the size of our data chunk actually is
      uint_t size    = cell.getData( function.getCellDataID() )->getSize( level );
      uint_t maxSize = getMaxDataBlockSize< Cell, func_t, value_t >( level );
      WALBERLA_CHECK_LESS_EQUAL( size, maxSize );

      // now we can set the position for this primitive's data chunk inside the global array
      WALBERLA_LOG_DETAIL( "Calling SetSelection for variable '" << varName << "' with\nstart .... {" << rowIndexInAdiosArray
                                                                 << ", 0}\n"
                                                                 << "count .... {1, " << size << "}" );
      varDoFData.SetSelection( { { rowIndexInAdiosArray, 0 }, { 1, size } } );

      // schedule the data for exporting
      engine.Put( varDoFData, cell.getData( function.getCellDataID() )->getPointer( level ) );

      break;
   }
   case Primitive::VERTEX:
      [[fallthrough]];
   case Primitive::EDGE:
      [[fallthrough]];
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type can only be FACE or CELL!" );
   }
}

template < template < typename > class func_t, concepts::value_type value_t >
   requires concepts::fe_function_scalar< func_t< value_t > >
void importVariableForScalarFunction( adios2::IO&              io,
                                      adios2::Engine&          engine,
                                      const func_t< value_t >& function,
                                      const std::string&       varName,
                                      uint_t                   level,
                                      uint_t                   rowIndexInAdiosArray,
                                      const Primitive&         primitive,
                                      uint_t                   step = 0U )
{
   WALBERLA_LOG_PROGRESS( "Running importVariableForScalarFunction() for '" << varName << "'" );

   // check that associated variable exists in checkpoint
   adios2::Variable< value_t > varDoFData = io.InquireVariable< value_t >( varName );
   if ( !varDoFData )
   {
      WALBERLA_ABORT( "ADIOS2 variable '" << varName << "' does not exist!" );
   }
   varDoFData.SetStepSelection( { step, 1u } );

   switch ( primitive.getType() )
   {
   case Primitive::FACE: {
      const Face& face = dynamic_cast< const Face& >( primitive );

      // figure out what the size of our data chunk actually is
      uint_t size    = face.getData( function.getFaceDataID() )->getSize( level );
      uint_t maxSize = getMaxDataBlockSize< Face, func_t, value_t >( level );
      WALBERLA_CHECK_LESS_EQUAL( size, maxSize );

      // now we can set the position for this primitive's data chunk inside the global array
      WALBERLA_LOG_DETAIL( "Calling SetSelection for variable '" << varName << "' with\nstart .... {" << rowIndexInAdiosArray
                                                                 << ", 0}\n"
                                                                 << "count .... {1, " << size << "}" );
      varDoFData.SetSelection( { { rowIndexInAdiosArray, 0 }, { 1, size } } );

      // schedule the data for import
      engine.Get( varDoFData, face.getData( function.getFaceDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::CELL: {
      const Cell& cell = dynamic_cast< const Cell& >( primitive );

      // figure out what the size of our data chunk actually is
      uint_t size    = cell.getData( function.getCellDataID() )->getSize( level );
      uint_t maxSize = getMaxDataBlockSize< Cell, func_t, value_t >( level );
      WALBERLA_CHECK_LESS_EQUAL( size, maxSize );

      // now we can set the position for this primitive's data chunk inside the global array
      WALBERLA_LOG_DETAIL( "Calling SetSelection for variable '" << varName << "' with\nstart .... {" << rowIndexInAdiosArray
                                                                 << ", 0}\n"
                                                                 << "count .... {1, " << size << "}" );
      varDoFData.SetSelection( { { rowIndexInAdiosArray, 0 }, { 1, size } } );

      engine.Get( varDoFData, cell.getData( function.getCellDataID() )->getPointer( level ) );
      break;
   }

   case Primitive::VERTEX:
      [[fallthrough]];
   case Primitive::EDGE:
      [[fallthrough]];
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type can only be FACE or CELL!" );
   }
}

template < template < typename > class func_t, concepts::value_type value_t >
void generateVariables( adios2::IO&              io,
                        const func_t< value_t >& function,
                        uint_t                   level,
                        uint_t                   globalNumberOfHighestDimPrimitives )
{
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                  std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         generateVariableForScalarFunction( io, function[k], level, globalNumberOfHighestDimPrimitives );
      }
   }

   else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                       std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         generateVariableForScalarFunction( io, function[k].getVertexDoFFunction(), level, globalNumberOfHighestDimPrimitives );
         generateVariableForScalarFunction( io, function[k].getEdgeDoFFunction(), level, globalNumberOfHighestDimPrimitives );
      }
   }

   else
   {
      WALBERLA_ABORT( "generateVariables() called with unsupported function type!" );
   }
}

template < template < typename > class func_t, concepts::value_type value_t >
   requires concepts::fe_function_scalar< func_t< value_t > > || concepts::fe_function_vectorial< func_t< value_t > >
void exportVariables( adios2::IO&              io,
                      adios2::Engine&          engine,
                      const func_t< value_t >& function,
                      uint_t                   level,
                      uint_t                   rowIndex,
                      const Primitive&         primitive,
                      uint_t                   step = 0U )
{
   // obtain row index in ADIOS2 array
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                  std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         exportVariableForScalarFunction( io, engine, function[k], level, rowIndex, primitive );
      }
   }

   else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                       std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         exportVariableForScalarFunction( io, engine, function[k].getVertexDoFFunction(), level, rowIndex, primitive );
         exportVariableForScalarFunction( io, engine, function[k].getEdgeDoFFunction(), level, rowIndex, primitive );
      }
   }

   else
   {
      WALBERLA_ABORT( "exportVariables() called with unsupported function type!" );
   }
}

template < template < typename > class func_t, concepts::value_type value_t >
   requires concepts::fe_function_scalar< func_t< value_t > > || concepts::fe_function_vectorial< func_t< value_t > >
void importVariables( adios2::IO&              io,
                      adios2::Engine&          engine,
                      const func_t< value_t >& function,
                      uint_t                   level,
                      uint_t                   rowIndex,
                      const Primitive&         primitive,
                      uint_t                   step = 0U )
{
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                  std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         std::string varName = generateVariableName( function[k].getFunctionName(), level );
         WALBERLA_LOG_PROGRESS( "Triggering import for '" << varName << "'" );
         importVariableForScalarFunction( io, engine, function[k], varName, level, rowIndex, primitive, step );
      }
   }

   else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                       std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         std::string varName = generateVariableName( function[k].getVertexDoFFunction().getFunctionName(), level );
         importVariableForScalarFunction(
             io, engine, function[k].getVertexDoFFunction(), varName, level, rowIndex, primitive, step );

         varName = generateVariableName( function[k].getEdgeDoFFunction().getFunctionName(), level );
         importVariableForScalarFunction(
             io, engine, function[k].getEdgeDoFFunction(), varName, level, rowIndex, primitive, step );
      }
   }
   else
   {
      WALBERLA_ABORT( "importVariables() called with unsupported function type!" );
   }
}

template < concepts::value_type value_t >
inline std::string valueTypeToString()
{
   std::string typeMarker;
   if constexpr ( std::is_same_v< value_t, double > )
   {
      typeMarker = "double";
   }
   else if constexpr ( std::is_same_v< value_t, float > )
   {
      typeMarker = "float";
   }
   else if constexpr ( std::is_same_v< value_t, int32_t > )
   {
      typeMarker = "int32_t";
   }
   else if constexpr ( std::is_same_v< value_t, int64_t > )
   {
      typeMarker = "int64_t";
   }
   else
   {
      WALBERLA_ABORT( "Achievement unlocked: 'Detector of the Missing Implementation'!" );
   }
   return typeMarker;
}

inline auto enumerateFaces( const std::shared_ptr< const PrimitiveStorage >& storage )
{
   // figure out some MPI runtime info
   uint_t myRank   = uint_c( walberla::MPIManager::instance()->rank() );
   uint_t commSize = uint_c( walberla::MPIManager::instance()->numProcesses() );

   // global communication to get face counts
   uint_t                numberOfLocalFaces = storage->getNumberOfLocalFaces();
   std::vector< uint_t > facesPerRank       = walberla::mpi::allGather( numberOfLocalFaces );

   // we enumerate w.r.t. MPI rank and local face counts
   uint_t offset = 0u;
   for ( uint_t i = 0; i < myRank; ++i )
   {
      offset += facesPerRank[i];
   }

   // enumerate local faces
   std::map< PrimitiveID, uint_t > localIndices;
   uint_t                          localCount = 0u;
   for ( const auto& item : storage->getFaces() )
   {
      localIndices[item.first] = offset + localCount;
      localCount++;
   }

   // figure out global number of faces
   uint_t globalNumberOfFaces = offset;
   for ( uint_t i = myRank; i < commSize; ++i )
   {
      globalNumberOfFaces += facesPerRank[i];
   }

   return std::make_pair( localIndices, globalNumberOfFaces );
}

inline auto enumerateCells( const std::shared_ptr< const PrimitiveStorage >& storage )
{
   // figure out some MPI runtime info
   uint_t myRank   = uint_c( walberla::MPIManager::instance()->rank() );
   uint_t commSize = uint_c( walberla::MPIManager::instance()->numProcesses() );

   // global communication to get cell counts
   uint_t                numberOfLocalCells = storage->getNumberOfLocalCells();
   std::vector< uint_t > cellsPerRank       = walberla::mpi::allGather( numberOfLocalCells );

   // we enumerate w.r.t. MPI rank and local face counts
   uint_t offset = 0u;
   for ( uint_t i = 0; i < myRank; ++i )
   {
      offset += cellsPerRank[i];
   }

   // enumerate local cells
   std::map< PrimitiveID, uint_t > localIndices;
   uint_t                          localCount = 0u;
   for ( const auto& item : storage->getCells() )
   {
      localIndices[item.first] = offset + localCount;
      localCount++;
   }

   // figure out global number of faces
   uint_t globalNumberOfCells = offset;
   for ( uint_t i = myRank; i < commSize; ++i )
   {
      globalNumberOfCells += cellsPerRank[i];
   }

   return std::make_pair( localIndices, globalNumberOfCells );
}

inline std::map< PrimitiveID, uint_t > deserializePrimitiveMap( const std::vector< uint8_t >& buffer, const adios2::Dims& shape )
{
   std::map< PrimitiveID, uint_t > primitiveMap;

   uint_t sizeOfKey   = sizeof( PrimitiveID );
   uint_t sizeOfValue = sizeof( uint_t );
   uint_t sizeOfPair  = sizeOfKey + sizeOfValue;

   // check that buffer size and number of pairs matches what we expect
   uint_t numPairs = buffer.size() / ( sizeOfKey + sizeOfValue );
   WALBERLA_CHECK( numPairs == shape[0] );
   WALBERLA_CHECK( ( sizeOfKey + sizeOfValue ) == shape[1] );

   const uint8_t* data = buffer.data();

   PrimitiveID key;
   uint_t      value;
   uint_t      bufferPos = 0u;
   for ( uint_t idx = 0; idx < numPairs; ++idx )
   {
      // extract key
      key.fromBuffer( buffer, bufferPos );
      bufferPos += sizeOfPair;

      // extract value
      value = *( reinterpret_cast< const uint_t* >( data + ( idx * sizeOfPair + sizeOfKey ) ) );

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Inserting pair into map:\n" << key << "\n" << value );

      // try insertion
      const auto [iterator, success] = primitiveMap.insert( { key, value } );
      if ( !success )
      {
         WALBERLA_ABORT( "Insertion of pair into map failed:\n" << key << "\n" << value );
      }
   }

   return primitiveMap;
}

template < concepts::value_type value_t >
void distributeVertexDoFData( const P1Function< value_t >& function, uint_t minLevel, uint_t maxLevel )
{
   if constexpr ( std::is_integral_v< value_t > )
   {
      WALBERLA_ABORT( "The current implementation of distributeVertexDoFData() does not support DoF data of integral type!" );
   }

   const auto& storage = function.getStorage();

   if ( storage->hasGlobalCells() )
   {
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         function.template communicateAdditively< Cell, Face >( level );
         function.template communicateAdditively< Cell, Edge >( level );
         function.template communicateAdditively< Cell, Vertex >( level );
      }
   }
   else
   {
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         function.template communicateAdditively< Face, Edge >( level, true );
         function.template communicateAdditively< Face, Vertex >( level, true );
      }
   }

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      if ( storage->hasGlobalCells() )
      {
         const auto& faceDataID = function.getFaceDataID();

         for ( const auto& item : storage->getFaces() )
         {
            Face& face = *( item.second );

            value_t factor = static_cast< value_t >( 1.0 / face.getNumNeighborCells() );
            vertexdof::macroface::assign( level, face, { factor }, { faceDataID }, faceDataID );
         }
      }

      const auto& edgeDataID   = function.getEdgeDataID();
      const auto& vertexDataID = function.getVertexDataID();

      for ( const auto& item : storage->getEdges() )
      {
         Edge& edge = *( item.second );

         value_t factor = static_cast< value_t >(
             1.0 / ( storage->hasGlobalCells() ? edge.getNumNeighborCells() : edge.getNumNeighborFaces() ) );
         vertexdof::macroedge::assign( level, edge, { factor }, { edgeDataID }, edgeDataID );
      }

      for ( const auto& item : storage->getVertices() )
      {
         Vertex& vertex = *( item.second );

         value_t factor = static_cast< value_t >(
             1.0 / ( storage->hasGlobalCells() ? vertex.getNumNeighborCells() : vertex.getNumNeighborFaces() ) );
         vertexdof::macrovertex::assign( vertex, { factor }, { vertexDataID }, vertexDataID, level );
      }
   }
}

template < concepts::value_type value_t >
void distributeEdgeDoFData( const EdgeDoFFunction< value_t >& function, uint_t minLevel, uint_t maxLevel )
{
   if constexpr ( std::is_integral_v< value_t > )
   {
      WALBERLA_ABORT( "The current implementation of distributeEdgeDoFData() does not support DoF data of integral type!" );
   }

   const auto& storage = function.getStorage();

   if ( storage->hasGlobalCells() )
   {
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         function.template communicateAdditively< Cell, Face >( level );
         function.template communicateAdditively< Cell, Edge >( level );
      }
   }
   else
   {
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         function.template communicateAdditively< Face, Edge >( level, true );
      }
   }

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      if ( storage->hasGlobalCells() )
      {
         const auto& faceDataID = function.getFaceDataID();

         for ( const auto& item : storage->getFaces() )
         {
            Face& face = *( item.second );

            value_t factor = static_cast< value_t >( 1.0 / face.getNumNeighborCells() );
            edgedof::macroface::assign( level, face, { factor }, { faceDataID }, faceDataID );
         }
      }

      const auto& edgeDataID = function.getEdgeDataID();

      for ( const auto& item : storage->getEdges() )
      {
         Edge& edge = *( item.second );

         value_t factor = static_cast< value_t >(
             1.0 / ( storage->hasGlobalCells() ? edge.getNumNeighborCells() : edge.getNumNeighborFaces() ) );
         edgedof::macroedge::assign( level, edge, { factor }, { edgeDataID }, edgeDataID );
      }
   }
}

template < template < typename > class func_t, concepts::value_type value_t >
   requires concepts::fe_function< func_t< value_t > >
void distributeDoFData( const func_t< value_t >& function, uint_t minLevel, uint_t maxLevel )
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Going to distribute DoF data to lower-dimensional primitives!" );

   const auto& storage = function.getStorage();

   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                  std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         distributeVertexDoFData( function[k], minLevel, maxLevel );
      }
   }

   else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                       std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         distributeVertexDoFData( function[k].getVertexDoFFunction(), minLevel, maxLevel );
         distributeEdgeDoFData( function[k].getEdgeDoFFunction(), minLevel, maxLevel );
      }
   }

   else
   {
      WALBERLA_ABORT( "importVariables() called with unsupported function type!" );
   }
}

} // namespace hyteg::adiosCheckpointHelpers_v03
