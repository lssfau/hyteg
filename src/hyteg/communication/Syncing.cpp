/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include "Syncing.hpp"

#include "hyteg/dataexport/FEFunctionRegistry.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
namespace communication {

using walberla::uint_t;

template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function, const uint_t& level )
{
   function.template communicate< Vertex, Edge >( level );
   function.template communicate< Edge, Face >( level );
   function.template communicate< Face, Cell >( level );

   function.template communicate< Cell, Face >( level );
   function.template communicate< Face, Edge >( level );
   function.template communicate< Edge, Vertex >( level );
}

// Version for VectorFunctions
template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< vType >& vecFunc, const uint_t& level )
{
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
   {
      syncFunctionBetweenPrimitives( vecFunc[idx], level );
   }
}

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< vType >& vecFunc, const uint_t& level )
{
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
   {
      syncFunctionBetweenPrimitives( vecFunc[idx], level );
   }
}

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const EGFunction< vType >& p1dgeFunc, const uint_t& level )
{
   auto& vecFunc = *( p1dgeFunc.getConformingPart() );
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
      syncFunctionBetweenPrimitives< hyteg::vertexdof::VertexDoFFunction< vType > >( vecFunc[idx], level );
}

template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< double >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< float >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< int32_t >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< int64_t >& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< double >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< float >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< int32_t >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< int64_t >& function, const uint_t& level );

template void syncFunctionBetweenPrimitives( const P2Function< double >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const P2Function< float >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const P2Function< int32_t >& function, const uint_t& level );
template void syncFunctionBetweenPrimitives( const P2Function< int64_t >& function, const uint_t& level );

template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< double >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< float >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< int32_t >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< int64_t >& function, const uint_t& level );

template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< double >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< float >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< int32_t >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< int64_t >& function, const uint_t& level );

template void syncVectorFunctionBetweenPrimitives( const EGFunction< double >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const EGFunction< int32_t >& function, const uint_t& level );
template void syncVectorFunctionBetweenPrimitives( const EGFunction< int64_t >& function, const uint_t& level );

void syncRegisteredFunctions( const FEFunctionRegistry& feFunctionRegistry, uint_t level )
{
   // ----------------------------------------
   //  P1Functions [double, int32_t, int64_t]
   // ----------------------------------------
   for ( const auto& function : feFunctionRegistry.getP1Functions().getFunctions< double >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives< hyteg::P1Function< double > >( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP1Functions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives< hyteg::P1Function< int32_t > >( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP1Functions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives< hyteg::P1Function< int64_t > >( function, level );
   }

   // ----------------------------------------------
   //  P1VectorFunctions [double, int32_t, int64_t]
   // ----------------------------------------------
   for ( const auto& function : feFunctionRegistry.getP1VectorFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP1VectorFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP1VectorFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }

   // ----------------------------------------
   //  P2Functions [double, int32_t, int64_t]
   // ----------------------------------------
   for ( const auto& function : feFunctionRegistry.getP2Functions().getFunctions< double >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP2Functions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP2Functions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level );
   }

   // ----------------------------------------------
   //  P2VectorFunctions [double, int32_t, int64_t]
   // ----------------------------------------------
   for ( const auto& function : feFunctionRegistry.getP2VectorFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP2VectorFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getP2VectorFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }

   // ----------------------------------------------
   //  DGVectorFunctions [double, int32_t, int64_t]
   // ----------------------------------------------

   // no communication necessary

   // ---------------------------------------------
   //  EdgeDoFFunctions [double, int32_t, int64_t]
   // ---------------------------------------------
   for ( const auto& function : feFunctionRegistry.getEdgeDoFFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getEdgeDoFFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getEdgeDoFFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level );
   }

   // ----------------------------------------
   //  DGFunctions [double, int32_t, int64_t]
   // ----------------------------------------

   // no communication necessary

   // -----------------------------------------------
   //  N1E1VectorFunction [double, int32_t, int64_t]
   // -----------------------------------------------
   for ( const auto& function : feFunctionRegistry.getN1E1VectorFunctions().getFunctions< double >() )
   {
      function.communicate< Face, Cell >( level );
      function.communicate< Edge, Cell >( level );
   }
   for ( const auto& function : feFunctionRegistry.getN1E1VectorFunctions().getFunctions< int32_t >() )
   {
      function.communicate< Face, Cell >( level );
      function.communicate< Edge, Cell >( level );
   }
   for ( const auto& function : feFunctionRegistry.getN1E1VectorFunctions().getFunctions< int64_t >() )
   {
      function.communicate< Face, Cell >( level );
      function.communicate< Edge, Cell >( level );
   }

   // ---------------------------------------------
   //  EGFunction [double, int32_t, int64_t]
   // ---------------------------------------------
   for ( const auto& function : feFunctionRegistry.getEGFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getEGFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }
   for ( const auto& function : feFunctionRegistry.getEGFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }
}

} // namespace communication
} // namespace hyteg
