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

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"

namespace hyteg::adiosCheckpointHelpers {

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
void generateVariableForScalarFunction( adios2::IO&              io,
                                        adios2::Engine&          engine,
                                        const func_t< value_t >& function,
                                        uint_t                   level,
                                        const Primitive&         primitive )
{
   WALBERLA_UNUSED( engine );
   WALBERLA_ASSERT( function.getDimension() == 1 );

   std::string varName = generateVariableName( function.getFunctionName(), primitive.getID(), level );
   uint_t      size    = 0;

   switch ( primitive.getType() )
   {
   case Primitive::VERTEX: {
      size = dynamic_cast< const Vertex& >( primitive ).getData( function.getVertexDataID() )->getSize( level );
      break;
   }
   case Primitive::EDGE: {
      size = dynamic_cast< const Edge& >( primitive ).getData( function.getEdgeDataID() )->getSize( level );
      break;
   }
   case Primitive::FACE: {
      size = dynamic_cast< const Face& >( primitive ).getData( function.getFaceDataID() )->getSize( level );
      break;
   }
   case Primitive::CELL: {
      size = dynamic_cast< const Cell& >( primitive ).getData( function.getCellDataID() )->getSize( level );
      break;
   }
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type is INVALID!" );
   }

   io.DefineVariable< value_t >( varName, {}, {}, { size } );
};

template < template < typename > class func_t, typename value_t >
void exportVariableForScalarFunction( adios2::IO&              io,
                                      adios2::Engine&          engine,
                                      const func_t< value_t >& function,
                                      uint_t                   level,
                                      const Primitive&         primitive )
{
   WALBERLA_ASSERT( function.getDimension() == 1 );

   std::string varName = generateVariableName( function.getFunctionName(), primitive.getID(), level );

   // check that associated variable exists in IO object
   adios2::Variable< value_t > varDoFData = io.InquireVariable< value_t >( varName );
   if ( !varDoFData )
   {
      WALBERLA_ABORT( "ADIOS2 variable '" << varName << "' does not exist!" );
   }

   switch ( primitive.getType() )
   {
   case Primitive::VERTEX: {
      const Vertex& vertex = dynamic_cast< const Vertex& >( primitive );
      engine.Put( varDoFData, vertex.getData( function.getVertexDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::EDGE: {
      const Edge& edge = dynamic_cast< const Edge& >( primitive );
      engine.Put( varDoFData, edge.getData( function.getEdgeDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::FACE: {
      const Face& face = dynamic_cast< const Face& >( primitive );
      engine.Put( varDoFData, face.getData( function.getFaceDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::CELL: {
      const Cell& cell = dynamic_cast< const Cell& >( primitive );
      engine.Put( varDoFData, cell.getData( function.getCellDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type is INVALID!" );
   }
}

template < template < typename > class func_t, typename value_t >
void importVariableForScalarFunction( adios2::IO&              io,
                                      adios2::Engine&          engine,
                                      const func_t< value_t >& function,
                                      const std::string&       varName,
                                      uint_t                   level,
                                      const Primitive&         primitive )
{
   WALBERLA_ASSERT( function.getDimension() == 1 );

   // check that associated variable exists in checkpoint
   adios2::Variable< value_t > varDoFData = io.InquireVariable< value_t >( varName );
   if ( !varDoFData )
   {
      WALBERLA_ABORT( "ADIOS2 variable '" << varName << "' does not exist!" );
   }

   switch ( primitive.getType() )
   {
   case Primitive::VERTEX: {
      const Vertex& vertex = dynamic_cast< const Vertex& >( primitive );
      engine.Get( varDoFData, vertex.getData( function.getVertexDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::EDGE: {
      const Edge& edge = dynamic_cast< const Edge& >( primitive );
      engine.Get( varDoFData, edge.getData( function.getEdgeDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::FACE: {
      const Face& face = dynamic_cast< const Face& >( primitive );
      engine.Get( varDoFData, face.getData( function.getFaceDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::CELL: {
      const Cell& cell = dynamic_cast< const Cell& >( primitive );
      engine.Get( varDoFData, cell.getData( function.getCellDataID() )->getPointer( level ) );
      break;
   }
   case Primitive::INVALID:
      WALBERLA_ABORT( "Primitive type is INVALID!" );
   }
}

template < template < typename > class func_t, typename value_t >
void generateVariables( adios2::IO&              io,
                        adios2::Engine&          engine,
                        const func_t< value_t >& function,
                        uint_t                   level,
                        const Primitive&         primitive )
{
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                  std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         generateVariableForScalarFunction( io, engine, function[k], level, primitive );
      }
   }

   else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                       std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         generateVariableForScalarFunction( io, engine, function[k].getVertexDoFFunction(), level, primitive );
         generateVariableForScalarFunction( io, engine, function[k].getEdgeDoFFunction(), level, primitive );
      }
   }

   else
   {
      WALBERLA_ABORT( "generateVariables() called with unsupported function type!" );
   }
}

template < template < typename > class func_t, typename value_t >
void exportVariables( adios2::IO&              io,
                      adios2::Engine&          engine,
                      const func_t< value_t >& function,
                      uint_t                   level,
                      const Primitive&         primitive )
{
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                  std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         exportVariableForScalarFunction( io, engine, function[k], level, primitive );
      }
   }

   else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                       std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         exportVariableForScalarFunction( io, engine, function[k].getVertexDoFFunction(), level, primitive );
         exportVariableForScalarFunction( io, engine, function[k].getEdgeDoFFunction(), level, primitive );
      }
   }

   else
   {
      WALBERLA_ABORT( "exportVariables() called with unsupported function type!" );
   }
}

template < template < typename > class func_t, typename value_t >
void importVariables( adios2::IO&              io,
                      adios2::Engine&          engine,
                      const func_t< value_t >& function,
                      uint_t                   level,
                      const Primitive&         primitive )
{
   if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                  std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         std::string varName = generateVariableName( function[k].getFunctionName(), primitive.getID(), level );
         importVariableForScalarFunction( io, engine, function[k], varName, level, primitive );
      }
   }

   else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                       std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
   {
      for ( uint_t k = 0; k < function.getDimension(); ++k )
      {
         std::string varName =
             generateVariableName( function[k].getVertexDoFFunction().getFunctionName(), primitive.getID(), level );
         importVariableForScalarFunction( io, engine, function[k].getVertexDoFFunction(), varName, level, primitive );

         varName = generateVariableName( function[k].getEdgeDoFFunction().getFunctionName(), primitive.getID(), level );
         importVariableForScalarFunction( io, engine, function[k].getEdgeDoFFunction(), varName, level, primitive );
      }
   }

   else
   {
      WALBERLA_ABORT( "importVariables() called with unsupported function type!" );
   }
}

template < typename value_t >
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

} // namespace hyteg::adiosCheckpointHelpers
