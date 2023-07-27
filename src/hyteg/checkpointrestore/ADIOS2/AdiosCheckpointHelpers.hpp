/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include "hyteg/PrimitiveID.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitives/Primitive.hpp"

namespace hyteg::adiosCheckpointHelpers {

template < typename value_t >
void exportDoFData( adios2::IO&                                io,
                    adios2::Engine&                            engine,
                    const std::string&                         varName,
                    const std::vector< const value_t const* >& dataBuffers,
                    const std::vector< uint_t >&               bufferSizes )
{
   WALBERLA_ASSERT( dataBuffers.size() == bufferSizes.size() );

   adios2::Variable< value_t > varDoFData = io.InquireVariable< value_t >( varName );
   if ( !varDoFData )
   {
      WALBERLA_ABORT( "ADIOS2 variable '" << varName << "' does not exist!" );
   }
   else
   {
      typename adios2::Variable< value_t >::Span dofData    = engine.Put< value_t >( varDoFData );
      uint_t                                     currentPos = 0;
      for ( uint_t k = 0; k < dataBuffers.size(); ++k )
      {
         memcpy( &dofData[currentPos], dataBuffers[k], bufferSizes[k] * sizeof( value_t ) );
         currentPos += bufferSizes[k];
      }
   }
};

std::string generateVariableName( const std::string& functionName, const PrimitiveID& id, uint_t level )
{
   std::stringstream sStream;
   sStream << functionName << "-ID";
   id.toStream( sStream );
   sStream << "-LVL" << level;
   return sStream.str();
}

template < template < typename > class func_t, typename value_t, typename action_t >
void doSomethingForAFunctionOnAllPrimitives( adios2::IO&              io,
                                             adios2::Engine&          engine,
                                             const func_t< value_t >& function,
                                             uint_t                   minLevel,
                                             uint_t                   maxLevel,
                                             action_t                 action )
{
   const auto& storage = function.getStorage();

   for ( const auto& macroIterator : storage->getCells() )
   {
      const Cell& cell = *macroIterator.second;
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         action( io, engine, function, level, cell );
      }
   }

   for ( const auto& macroIterator : storage->getFaces() )
   {
      const Face& face = *macroIterator.second;
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         action( io, engine, function, level, face );
      }
   }

   for ( const auto& macroIterator : storage->getEdges() )
   {
      const Edge& edge = *macroIterator.second;
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         action( io, engine, function, level, edge );
      }
   }

   for ( const auto& macroIterator : storage->getVertices() )
   {
      const Vertex& vertex = *macroIterator.second;
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         action( io, engine, function, level, vertex );
      }
   }
}

template < template < typename > class func_t, typename value_t >
void generateVariableForP1TypeFunction( adios2::IO&              io,
                                        adios2::Engine&          engine,
                                        const func_t< value_t >& function,
                                        uint_t                   level,
                                        const Primitive&         primitive )
{
   WALBERLA_UNUSED( engine );

   if constexpr ( !std::is_same_v< func_t< value_t >, P1Function< value_t > > &&
                  !std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      WALBERLA_ABORT( "generateVariableForP1TypeFunction() called with a non P1 kind function!" );
   }

   std::string varName = generateVariableName( function.getFunctionName(), primitive.getID(), level );

   // We cannot ask a P1VectorFunction for a DataID (needed below), but only its component functions,
   // so we do a little abstraction here, to minimise code
   const P1Function< value_t > const* functionPtr = nullptr;
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > )
   {
      functionPtr = std::addressof( function );
   }
   else
   {
      functionPtr = std::addressof( function[0] );
   }

   switch ( primitive.getType() )
   {
   case Primitive::VERTEX: {
      uint_t size = dynamic_cast< const Vertex& >( primitive ).getData( functionPtr->getVertexDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::EDGE: {
      uint_t size = dynamic_cast< const Edge& >( primitive ).getData( functionPtr->getEdgeDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::FACE: {
      uint_t size = dynamic_cast< const Face& >( primitive ).getData( functionPtr->getFaceDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::CELL: {
      uint_t size = dynamic_cast< const Cell& >( primitive ).getData( functionPtr->getCellDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type is INVALID!" );
   }

   // WALBERLA_LOG_INFO_ON_ROOT( "Defined variable '" << varName << "'" );
};

template < template < typename > class func_t, typename value_t >
void exportVariableForP1TypeFunction( adios2::IO&              io,
                                      adios2::Engine&          engine,
                                      const func_t< value_t >& function,
                                      uint_t                   level,
                                      const Primitive&         primitive )
{
   if constexpr ( !std::is_same_v< func_t< value_t >, P1Function< value_t > > &&
                  !std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      WALBERLA_ABORT( "exportVariableForP1TypeFunction() called with a non-P1-type function!" );
   }

   // check that associated variable exists in IO object
   std::string                 varName    = generateVariableName( function.getFunctionName(), primitive.getID(), level );
   adios2::Variable< value_t > varDoFData = io.InquireVariable< value_t >( varName );
   if ( !varDoFData )
   {
      WALBERLA_ABORT( "ADIOS2 variable '" << varName << "' does not exist!" );
   }

   // For a P1VectorFunction we need to work on its component functions, so we do a little abstraction here,
   // to minimise coding
   std::vector< const P1Function< value_t > const* > componentFunction;
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > )
   {
      componentFunction.push_back( std::addressof( function ) );
   }
   else
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         componentFunction.push_back( std::addressof( function[k] ) );
      }
   }

   // fetch span to copy our data into
   typename adios2::Variable< value_t >::Span dofData = engine.Put< value_t >( varDoFData );

   switch ( primitive.getType() )
   {
   case Primitive::VERTEX: {
      const Vertex& vertex     = dynamic_cast< const Vertex& >( primitive );
      uint_t        size       = vertex.getData( componentFunction[0]->getVertexDataID() )->getSize( level );
      uint_t        currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer = vertex.getData( componentFunction[k]->getVertexDataID() )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer, size * sizeof( value_t ) );
         currentPos += size;
      }
      break;
   }
   case Primitive::EDGE: {
      const Edge& edge       = dynamic_cast< const Edge& >( primitive );
      uint_t      size       = edge.getData( componentFunction[0]->getEdgeDataID() )->getSize( level );
      uint_t      currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer = edge.getData( componentFunction[k]->getEdgeDataID() )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer, size * sizeof( value_t ) );
         currentPos += size;
      }
      break;
   }
   case Primitive::FACE: {
      const Face& face       = dynamic_cast< const Face& >( primitive );
      uint_t      size       = face.getData( componentFunction[0]->getFaceDataID() )->getSize( level );
      uint_t      currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer = face.getData( componentFunction[k]->getFaceDataID() )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer, size * sizeof( value_t ) );
         currentPos += size;
      }
      break;
   }
   case Primitive::CELL: {
      const Cell& cell       = dynamic_cast< const Cell& >( primitive );
      uint_t      size       = cell.getData( componentFunction[0]->getCellDataID() )->getSize( level );
      uint_t      currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer = cell.getData( componentFunction[k]->getCellDataID() )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer, size * sizeof( value_t ) );
         currentPos += size;
      }
      break;
   }
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type is INVALID!" );
   }
};

template < template < typename > class func_t, typename value_t >
void generateVariableForP2TypeFunction( adios2::IO&              io,
                                        adios2::Engine&          engine,
                                        const func_t< value_t >& function,
                                        uint_t                   level,
                                        const Primitive&         primitive )
{
   WALBERLA_UNUSED( engine );

   if constexpr ( !std::is_same_v< func_t< value_t >, P2Function< value_t > > &&
                  !std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      WALBERLA_ABORT( "generateVariableForP2TypeFunction() called with a non-P2-type function!" );
   }

   std::string varName = generateVariableName( function.getFunctionName(), primitive.getID(), level );

   // We cannot ask a P2VectorFunction for a DataID (needed below), but only its component functions,
   // so we do a little abstraction here, to minimise code
   const P2Function< value_t > const* functionPtr = nullptr;
   if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > )
   {
      functionPtr = std::addressof( function );
   }
   else
   {
      functionPtr = std::addressof( function[0] );
   }

   switch ( primitive.getType() )
   {
   case Primitive::VERTEX: {
      const Vertex& vertex = dynamic_cast< const Vertex& >( primitive );
      uint_t        size   = vertex.getData( functionPtr->getVertexDoFFunction().getVertexDataID() )->getSize( level );
      size += vertex.getData( functionPtr->getEdgeDoFFunction().getVertexDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::EDGE: {
      const Edge& edge = dynamic_cast< const Edge& >( primitive );
      uint_t      size = edge.getData( functionPtr->getVertexDoFFunction().getEdgeDataID() )->getSize( level );
      size += edge.getData( functionPtr->getEdgeDoFFunction().getEdgeDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::FACE: {
      const Face& face = dynamic_cast< const Face& >( primitive );
      uint_t      size = face.getData( functionPtr->getVertexDoFFunction().getFaceDataID() )->getSize( level );
      size += face.getData( functionPtr->getEdgeDoFFunction().getFaceDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::CELL: {
      const Cell& cell = dynamic_cast< const Cell& >( primitive );
      uint_t      size = cell.getData( functionPtr->getVertexDoFFunction().getCellDataID() )->getSize( level );
      size += cell.getData( functionPtr->getEdgeDoFFunction().getCellDataID() )->getSize( level );
      io.DefineVariable< value_t >( varName, {}, {}, { size * function.getDimension() } );
      break;
   }
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type is INVALID!" );
   }

   // WALBERLA_LOG_INFO_ON_ROOT( "Defined variable '" << varName << "'" );
};

template < template < typename > class func_t, typename value_t >
void exportVariableForP2TypeFunction( adios2::IO&              io,
                                      adios2::Engine&          engine,
                                      const func_t< value_t >& function,
                                      uint_t                   level,
                                      const Primitive&         primitive )
{
   if constexpr ( !std::is_same_v< func_t< value_t >, P2Function< value_t > > &&
                  !std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      WALBERLA_ABORT( "exportVariableForP2TypeFunction() called with a non-P2-type function!" );
   }

   // check that associated variable exists in IO object
   std::string                 varName    = generateVariableName( function.getFunctionName(), primitive.getID(), level );
   adios2::Variable< value_t > varDoFData = io.InquireVariable< value_t >( varName );
   if ( !varDoFData )
   {
      WALBERLA_ABORT( "ADIOS2 variable  '" << varName << "' does not exist!" );
   }

   // For a P2VectorFunction we need to work on its component functions, so we do a little abstraction here,
   // to minimise coding
   std::vector< const P2Function< value_t > const* > componentFunction;
   if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > )
   {
      componentFunction.push_back( std::addressof( function ) );
   }
   else
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         componentFunction.push_back( std::addressof( function[k] ) );
      }
   }

   // fetch span to copy our data into
   typename adios2::Variable< value_t >::Span dofData = engine.Put< value_t >( varDoFData );

   switch ( primitive.getType() )
   {
   case Primitive::VERTEX: {
      const Vertex& vertex      = dynamic_cast< const Vertex& >( primitive );
      const auto    VertexDofID = componentFunction[0]->getVertexDoFFunction().getVertexDataID();
      const auto    EdgeDofID   = componentFunction[0]->getEdgeDoFFunction().getVertexDataID();

      uint_t size1      = vertex.getData( VertexDofID )->getSize( level );
      uint_t size2      = vertex.getData( EdgeDofID )->getSize( level );
      uint_t currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer1 = vertex.getData( VertexDofID )->getPointer( level );
         value_t* dataBuffer2 = vertex.getData( EdgeDofID )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer1, size1 * sizeof( value_t ) );
         currentPos += size1;
         memcpy( &dofData[currentPos], dataBuffer2, size2 * sizeof( value_t ) );
         currentPos += size2;
      }
      break;
   }
   case Primitive::EDGE: {
      const Edge& edge        = dynamic_cast< const Edge& >( primitive );
      const auto  VertexDofID = componentFunction[0]->getVertexDoFFunction().getEdgeDataID();
      const auto  EdgeDofID   = componentFunction[0]->getEdgeDoFFunction().getEdgeDataID();

      uint_t size1      = edge.getData( VertexDofID )->getSize( level );
      uint_t size2      = edge.getData( EdgeDofID )->getSize( level );
      uint_t currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer1 = edge.getData( VertexDofID )->getPointer( level );
         value_t* dataBuffer2 = edge.getData( EdgeDofID )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer1, size1 * sizeof( value_t ) );
         currentPos += size1;
         memcpy( &dofData[currentPos], dataBuffer2, size2 * sizeof( value_t ) );
         currentPos += size2;
      }
      break;
   }
   case Primitive::FACE: {
      const Face& face        = dynamic_cast< const Face& >( primitive );
      const auto  VertexDofID = componentFunction[0]->getVertexDoFFunction().getFaceDataID();
      const auto  EdgeDofID   = componentFunction[0]->getEdgeDoFFunction().getFaceDataID();

      uint_t size1      = face.getData( VertexDofID )->getSize( level );
      uint_t size2      = face.getData( EdgeDofID )->getSize( level );
      uint_t currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer1 = face.getData( VertexDofID )->getPointer( level );
         value_t* dataBuffer2 = face.getData( EdgeDofID )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer1, size1 * sizeof( value_t ) );
         currentPos += size1;
         memcpy( &dofData[currentPos], dataBuffer2, size2 * sizeof( value_t ) );
         currentPos += size2;
      }
      break;
   }
   case Primitive::CELL: {
      const Cell& cell        = dynamic_cast< const Cell& >( primitive );
      const auto  VertexDofID = componentFunction[0]->getVertexDoFFunction().getCellDataID();
      const auto  EdgeDofID   = componentFunction[0]->getEdgeDoFFunction().getCellDataID();

      uint_t size1      = cell.getData( VertexDofID )->getSize( level );
      uint_t size2      = cell.getData( EdgeDofID )->getSize( level );
      uint_t currentPos = 0;

      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         value_t* dataBuffer1 = cell.getData( VertexDofID )->getPointer( level );
         value_t* dataBuffer2 = cell.getData( EdgeDofID )->getPointer( level );
         memcpy( &dofData[currentPos], dataBuffer1, size1 * sizeof( value_t ) );
         currentPos += size1;
         memcpy( &dofData[currentPos], dataBuffer2, size2 * sizeof( value_t ) );
         currentPos += size2;
      }
      break;
   }
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type is INVALID!" );
   }
};

} // namespace hyteg::adiosCheckpointHelpers
