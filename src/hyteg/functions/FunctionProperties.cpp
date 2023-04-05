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

#include "hyteg/functions/FunctionProperties.hpp"

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"

namespace hyteg {

void printFunctionAllocationInfo( const PrimitiveStorage& storage, const uint_t& verbosityLevel )
{
   const uint_t sizeRealT = sizeof( real_t );

   const auto vertexDoFFunctionLevelWiseCounter =
       Function< vertexdof::VertexDoFFunction< real_t > >::getLevelWiseFunctionCounter();
   const auto edgeDoFFunctionLevelWiseCounter = Function< EdgeDoFFunction< real_t > >::getLevelWiseFunctionCounter();

   // number of functions
   const uint_t numVertexDoFFunctions = Function< vertexdof::VertexDoFFunction< real_t > >::getNumFunctions();
   const uint_t numEdgeDoFFunctions   = Function< EdgeDoFFunction< real_t > >::getNumFunctions();

   // get min and max level globally
   uint_t minLevelVertexDoF = std::numeric_limits< uint_t >::max();
   uint_t minLevelEdgeDoF   = std::numeric_limits< uint_t >::max();
   uint_t maxLevelVertexDoF = 0;
   uint_t maxLevelEdgeDoF   = 0;

   for ( auto it : vertexDoFFunctionLevelWiseCounter )
   {
      if ( it.first < minLevelVertexDoF )
         minLevelVertexDoF = it.first;
      if ( it.first > maxLevelVertexDoF )
         maxLevelVertexDoF = it.first;
   }

   for ( auto it : edgeDoFFunctionLevelWiseCounter )
   {
      if ( it.first < minLevelEdgeDoF )
         minLevelEdgeDoF = it.first;
      if ( it.first > maxLevelEdgeDoF )
         maxLevelEdgeDoF = it.first;
   }

   // number of DoFs
   uint_t numVertexDoFs = 0;
   uint_t numEdgeDoFs   = 0;
   for ( auto it : vertexDoFFunctionLevelWiseCounter )
   {
      const auto level        = it.first;
      const auto numFunctions = it.second;
      numVertexDoFs += numFunctions * numberOfGlobalDoFs< VertexDoFFunctionTag >( storage, level );
   }
   for ( auto it : edgeDoFFunctionLevelWiseCounter )
   {
      const auto level        = it.first;
      const auto numFunctions = it.second;
      numEdgeDoFs += numFunctions * numberOfGlobalDoFs< EdgeDoFFunctionTag >( storage, level );
   }

   // memory per DoF
   const uint_t vertexdofBytes = sizeRealT * numVertexDoFs;
   const double vertexdofGB    = real_c( vertexdofBytes ) / 1e+09;
   const uint_t edgedofBytes   = sizeRealT * numEdgeDoFs;
   const double edgedofGB      = real_c( edgedofBytes ) / 1e+09;

   // calculate actual allocated memory
   const double globalActualAllocatedMemory = double( FunctionMemory< real_t >::getGlobalAllocatedMemoryInBytes() ) / 1e+09;
   const double minActualAllocatedMemory    = double( FunctionMemory< real_t >::getMinLocalAllocatedMemoryInBytes() ) / 1e+09;
   const double maxActualAllocatedMemory    = double( FunctionMemory< real_t >::getMaxLocalAllocatedMemoryInBytes() ) / 1e+09;

   if ( verbosityLevel == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "[Function] allocated memory (sum | min | max) in GB: "
                                 << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << globalActualAllocatedMemory
                                 << " | " << std::fixed << std::setw( 12 ) << std::setprecision( 3 ) << minActualAllocatedMemory
                                 << " | " << std::fixed << std::setw( 12 ) << std::setprecision( 3 ) << maxActualAllocatedMemory )
   }

   if ( verbosityLevel > 0 )
   {
      // print it all...
      WALBERLA_LOG_INFO_ON_ROOT( "====================== Function Allocation Info ======================" );
      WALBERLA_LOG_INFO_ON_ROOT( "                       +--------------+--------------+--------------+" );
      WALBERLA_LOG_INFO_ON_ROOT( "                       |          sum |          min |          max |" );
      WALBERLA_LOG_INFO_ON_ROOT( " ----------------------+--------------+--------------+--------------+" );
      WALBERLA_LOG_INFO_ON_ROOT( " allocated memory (GB) | "
                                 << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << globalActualAllocatedMemory
                                 << " | " << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << minActualAllocatedMemory
                                 << " | " << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << maxActualAllocatedMemory
                                 << " | " );
      WALBERLA_LOG_INFO_ON_ROOT( " ----------------------+--------------+--------------+--------------+" );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      WALBERLA_LOG_INFO_ON_ROOT( "           +----------------------+----------------------+----------------------+" );
      WALBERLA_LOG_INFO_ON_ROOT( "           |  number of functions |          number DoFs | memory DoFs only (GB)|" );
      WALBERLA_LOG_INFO_ON_ROOT( " ----------+----------------------+----------------------+----------------------+" );
      WALBERLA_LOG_INFO_ON_ROOT( " VertexDoF | " << std::setw( 20 ) << numVertexDoFFunctions << " | " << std::setw( 20 )
                                                 << numVertexDoFs << " | " << std::setw( 20 ) << std::fixed
                                                 << std::setprecision( 3 ) << vertexdofGB << " |" );
      WALBERLA_LOG_INFO_ON_ROOT( " EdgeDoF   | " << std::setw( 20 ) << numEdgeDoFFunctions << " | " << std::setw( 20 )
                                                 << numEdgeDoFs << " | " << std::setw( 20 ) << std::fixed
                                                 << std::setprecision( 3 ) << edgedofGB << " |" );
      WALBERLA_LOG_INFO_ON_ROOT( " ----------+----------------------+----------------------+----------------------+" );
      WALBERLA_LOG_INFO_ON_ROOT( " total     | " << std::setw( 20 ) << numVertexDoFFunctions + numEdgeDoFFunctions << " | "
                                                 << std::setw( 20 ) << numVertexDoFs + numEdgeDoFs << " | " << std::setw( 20 )
                                                 << std::fixed << std::setprecision( 3 ) << vertexdofGB + edgedofGB << " |" );
      WALBERLA_LOG_INFO_ON_ROOT( " ----------+----------------------+----------------------+----------------------+" );

      WALBERLA_LOG_INFO_ON_ROOT( "" )
      WALBERLA_LOG_INFO_ON_ROOT( " - mixed and composites:" );
      WALBERLA_LOG_INFO_ON_ROOT( "   + P2Function:       " << Function< P2Function< real_t > >::getNumFunctions() );

      if ( verbosityLevel > 1 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "" )
         WALBERLA_LOG_INFO_ON_ROOT( " - currently allocated functions:" )
         WALBERLA_LOG_INFO_ON_ROOT( "   + VertexDoF (" << numVertexDoFFunctions << "):" )
         for ( auto name : Function< vertexdof::VertexDoFFunction< real_t > >::getFunctionNames() )
            WALBERLA_LOG_INFO_ON_ROOT( "      - " << name )
         WALBERLA_LOG_INFO_ON_ROOT( "   + EdgeDoF (" << numEdgeDoFFunctions << "):" )
         for ( auto name : Function< EdgeDoFFunction< real_t > >::getFunctionNames() )
            WALBERLA_LOG_INFO_ON_ROOT( "      - " << name )
         WALBERLA_LOG_INFO_ON_ROOT( "   + P2 (" << Function< P2Function< real_t > >::getNumFunctions() << "):" )
         for ( auto name : Function< P2Function< real_t > >::getFunctionNames() )
            WALBERLA_LOG_INFO_ON_ROOT( "      - " << name )
      }

      WALBERLA_LOG_INFO_ON_ROOT( "======================================================================" );
   }
}

} // namespace hyteg
