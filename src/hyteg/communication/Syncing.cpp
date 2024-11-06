/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/functions/FEFunctionRegistry.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
namespace communication {

using walberla::uint_t;

template < typename funcType >
void syncFunctionBetweenPrimitives( const funcType& function, const uint_t& level, syncDirection_t direction )
{
   if ( direction == syncDirection_t::LOW2HIGH )
   {
      function.template communicate< Vertex, Edge >( level );
      function.template communicate< Edge, Face >( level );
      function.template communicate< Face, Cell >( level );
   }

   else if ( direction == syncDirection_t::BIDIRECTIONAL )
   {
      function.template communicate< Vertex, Edge >( level );
      function.template communicate< Edge, Face >( level );
      function.template communicate< Face, Cell >( level );

      function.template communicate< Cell, Face >( level );
      function.template communicate< Face, Edge >( level );
      function.template communicate< Edge, Vertex >( level );
   }
}

// Version for VectorFunctions
template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< vType >& vecFunc,
                                          const uint_t&                    level,
                                          syncDirection_t                  direction )
{
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
   {
      syncFunctionBetweenPrimitives( vecFunc[idx], level, direction );
   }
}

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< vType >& vecFunc,
                                          const uint_t&                    level,
                                          syncDirection_t                  direction )
{
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
   {
      syncFunctionBetweenPrimitives( vecFunc[idx], level, direction );
   }
}

template < typename vType >
void syncVectorFunctionBetweenPrimitives( const EGFunction< vType >& p1dgeFunc, const uint_t& level, syncDirection_t direction )
{
   auto& vecFunc = *( p1dgeFunc.getConformingPart() );
   for ( uint_t idx = 0; idx < vecFunc.getDimension(); ++idx )
      syncFunctionBetweenPrimitives< hyteg::vertexdof::VertexDoFFunction< vType > >( vecFunc[idx], level, direction );
}

// -----------------------------------------------
//  Explicit Instantiations for VertexDoFFunction
// -----------------------------------------------
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< double >& function,
                                             const uint_t&                                 level,
                                             syncDirection_t                               direction );

template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< float >& function,
                                             const uint_t&                                level,
                                             syncDirection_t                              direction );

#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< walberla::float16 >& function,
                                             const uint_t&                                            level,
                                             syncDirection_t                                          direction );
#endif

template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< int32_t >& function,
                                             const uint_t&                                  level,
                                             syncDirection_t                                direction );

template void syncFunctionBetweenPrimitives( const vertexdof::VertexDoFFunction< int64_t >& function,
                                             const uint_t&                                  level,
                                             syncDirection_t                                direction );

// ---------------------------------------------
//  Explicit Instantiations for EdgeDoFFunction
// ---------------------------------------------
template void
    syncFunctionBetweenPrimitives( const EdgeDoFFunction< double >& function, const uint_t& level, syncDirection_t direction );

template void
    syncFunctionBetweenPrimitives( const EdgeDoFFunction< float >& function, const uint_t& level, syncDirection_t direction );

#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void syncFunctionBetweenPrimitives( const EdgeDoFFunction< walberla::float16 >& function,
                                             const uint_t&                               level,
                                             syncDirection_t                             direction );
#endif

template void
    syncFunctionBetweenPrimitives( const EdgeDoFFunction< int32_t >& function, const uint_t& level, syncDirection_t direction );

template void
    syncFunctionBetweenPrimitives( const EdgeDoFFunction< int64_t >& function, const uint_t& level, syncDirection_t direction );

// ----------------------------------------
//  Explicit Instantiations for P2Function
// ----------------------------------------
template void
    syncFunctionBetweenPrimitives( const P2Function< double >& function, const uint_t& level, syncDirection_t direction );

template void
    syncFunctionBetweenPrimitives( const P2Function< float >& function, const uint_t& level, syncDirection_t direction );

template void
    syncFunctionBetweenPrimitives( const P2Function< int32_t >& function, const uint_t& level, syncDirection_t direction );

template void
    syncFunctionBetweenPrimitives( const P2Function< int64_t >& function, const uint_t& level, syncDirection_t direction );

// --------------------------------------------------
//  Explicit Instantiations for P2PlusBubbleFunction
// --------------------------------------------------
template void syncFunctionBetweenPrimitives( const P2PlusBubbleFunction< double >& function,
                                             const uint_t&                         level,
                                             syncDirection_t                       direction );

template void syncFunctionBetweenPrimitives( const P2PlusBubbleFunction< float >& function,
                                             const uint_t&                        level,
                                             syncDirection_t                      direction );

template void syncFunctionBetweenPrimitives( const P2PlusBubbleFunction< int32_t >& function,
                                             const uint_t&                          level,
                                             syncDirection_t                        direction );

template void syncFunctionBetweenPrimitives( const P2PlusBubbleFunction< int64_t >& function,
                                             const uint_t&                          level,
                                             syncDirection_t                        direction );

// ----------------------------------------------
//  Explicit Instantiations for P1VectorFunction
// ----------------------------------------------
template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< double >& function,
                                                   const uint_t&                     level,
                                                   syncDirection_t                   direction );

template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< float >& function,
                                                   const uint_t&                    level,
                                                   syncDirection_t                  direction );

template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< int32_t >& function,
                                                   const uint_t&                      level,
                                                   syncDirection_t                    direction );

template void syncVectorFunctionBetweenPrimitives( const P1VectorFunction< int64_t >& function,
                                                   const uint_t&                      level,
                                                   syncDirection_t                    direction );

// ----------------------------------------------
//  Explicit Instantiations for P2VectorFunction
// ----------------------------------------------
template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< double >& function,
                                                   const uint_t&                     level,
                                                   syncDirection_t                   direction );

template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< float >& function,
                                                   const uint_t&                    level,
                                                   syncDirection_t                  direction );

template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< int32_t >& function,
                                                   const uint_t&                      level,
                                                   syncDirection_t                    direction );

template void syncVectorFunctionBetweenPrimitives( const P2VectorFunction< int64_t >& function,
                                                   const uint_t&                      level,
                                                   syncDirection_t                    direction );

// ----------------------------------------
//  Explicit Instantiations for EGFunction
// ----------------------------------------
template void
    syncVectorFunctionBetweenPrimitives( const EGFunction< double >& function, const uint_t& level, syncDirection_t direction );

template void
    syncVectorFunctionBetweenPrimitives( const EGFunction< float >& function, const uint_t& level, syncDirection_t direction );

template void
    syncVectorFunctionBetweenPrimitives( const EGFunction< int32_t >& function, const uint_t& level, syncDirection_t direction );

template void
    syncVectorFunctionBetweenPrimitives( const EGFunction< int64_t >& function, const uint_t& level, syncDirection_t direction );

void syncRegisteredFunctions( const FEFunctionRegistry& feFunctionRegistry,
                              uint_t                    level,
                              bool                      excludeDGTypeFunctions,
                              syncDirection_t           direction )
{
   uint_t controlCount{ 0u };

   // -----------------------------------------------
   //  P1Functions [double, float, int32_t, int64_t]
   // -----------------------------------------------
   for ( const auto& function : feFunctionRegistry.getP1Functions().getFunctions< double >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP1Functions().getFunctions< float >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP1Functions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP1Functions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }

   // -----------------------------------------------------
   //  P1VectorFunctions [double, float, int32_t, int64_t]
   // -----------------------------------------------------
   for ( const auto& function : feFunctionRegistry.getP1VectorFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP1VectorFunctions().getFunctions< float >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP1VectorFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP1VectorFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }

   // -----------------------------------------------
   //  P2Functions [double, float, int32_t, int64_t]
   // -----------------------------------------------
   for ( const auto& function : feFunctionRegistry.getP2Functions().getFunctions< double >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2Functions().getFunctions< float >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2Functions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2Functions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }

   // -----------------------------------------------------
   //  P2VectorFunctions [double, float, int32_t, int64_t]
   // -----------------------------------------------------
   for ( const auto& function : feFunctionRegistry.getP2VectorFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2VectorFunctions().getFunctions< float >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2VectorFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2VectorFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }

   // ---------------------------------------------------------
   //  P2PlusBubbleFunctions [double, float, int32_t, int64_t]
   // ---------------------------------------------------------
   if ( !excludeDGTypeFunctions && feFunctionRegistry.getP2PlusBubbleFunctions().size() > 0 )
   {
      WALBERLA_ABORT( "Sorry, but P2PlusBubbleFunction::communicate() currently does not communicate bubble dofs!\n"
                      << "Make sure to call syncing with excludeDGTypeFunctions = true!" );
   }
   for ( const auto& function : feFunctionRegistry.getP2PlusBubbleFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2PlusBubbleFunctions().getFunctions< float >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2PlusBubbleFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getP2PlusBubbleFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }

   // ---------------------------------------------
   //  DGFunctions [double, foat,int32_t, int64_t]
   // ---------------------------------------------
   if ( !excludeDGTypeFunctions && feFunctionRegistry.getDGFunctions().size() > 0 )
   {
      WALBERLA_ABORT( "Sorry implementation missing in syncRegisteredFunctions() for DGFunction objects!" );
   }
   controlCount += feFunctionRegistry.getDGFunctions().size();

   // ---------------------------------------------------
   //  DGVectorFunctions [double, foat,int32_t, int64_t]
   // ---------------------------------------------------
   if ( !excludeDGTypeFunctions )
   {
      WALBERLA_ABORT( "Sorry implementation missing in syncRegisteredFunctions() for DGVectorFunction objects!" );
   }
   controlCount += feFunctionRegistry.getDGVectorFunctions().size();

   // ----------------------------------------------------
   //  EdgeDoFFunctions [double, float, int32_t, int64_t]
   // ----------------------------------------------------
   for ( const auto& function : feFunctionRegistry.getEdgeDoFFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getEdgeDoFFunctions().getFunctions< float >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getEdgeDoFFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getEdgeDoFFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }

   // ------------------------------------------------------
   //  N1E1VectorFunction [double, float, int32_t, int64_t]
   // ------------------------------------------------------
   for ( const auto& function : feFunctionRegistry.getN1E1VectorFunctions().getFunctions< double >() )
   {
      function.communicate< Face, Cell >( level );
      function.communicate< Edge, Cell >( level );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getN1E1VectorFunctions().getFunctions< float >() )
   {
      function.communicate< Face, Cell >( level );
      function.communicate< Edge, Cell >( level );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getN1E1VectorFunctions().getFunctions< int32_t >() )
   {
      function.communicate< Face, Cell >( level );
      function.communicate< Edge, Cell >( level );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getN1E1VectorFunctions().getFunctions< int64_t >() )
   {
      function.communicate< Face, Cell >( level );
      function.communicate< Edge, Cell >( level );
      controlCount++;
   }

   // ----------------------------------------------
   //  EGFunction [double, float, int32_t, int64_t]
   // ----------------------------------------------
   for ( const auto& function : feFunctionRegistry.getEGFunctions().getFunctions< double >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getEGFunctions().getFunctions< float >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getEGFunctions().getFunctions< int32_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }
   for ( const auto& function : feFunctionRegistry.getEGFunctions().getFunctions< int64_t >() )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level, direction );
      controlCount++;
   }

   WALBERLA_CHECK_EQUAL( controlCount,
                         feFunctionRegistry.getNumberOfFunctionsInRegistry(),
                         "Synced " << controlCount << " functions, but registry stores "
                                   << feFunctionRegistry.getNumberOfFunctionsInRegistry()
                                   << ". Maybe there is a valueType template instance missing?" );
}

} // namespace communication
} // namespace hyteg
