/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionIterator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

void testWeightsInCell( const uint_t& lowerLevel )
{
   typedef edgedof::EdgeDoFOrientation eo;

   const auto storage = PrimitiveStorage::createFromGmshFile( prependHyTeGMeshDir( "3D/tet_1el.msh" ) );

   P2Function< real_t > u( "u", storage, lowerLevel, lowerLevel + 1 );
   u.interpolate( 1.0, lowerLevel + 1 );

   P2toP2QuadraticRestriction restrictionOperator;
   restrictionOperator.restrict( u, lowerLevel + 1, All );

   for ( auto it : FunctionIterator< P1Function< real_t > >( u.getVertexDoFFunction(), lowerLevel ) )
   {
      if ( it.isOnMacroCell() )
      {
         if ( it.isVertexDoF() )
         {
            // 24 neighbor cells
            // 14 neighbor edges
            real_t expected = 0;
            expected += 1.0;                         // vertex
            expected += 14.0 * 3.0 / 8.0;            // near edgedofs
            expected += 14.0 * ( -1.0 / 8.0 );       // far edgedofs
            expected += 3.0 * 24.0 * ( -1.0 / 8.0 ); // edgedofs on faces of cell
            expected += 24.0 * ( -1.0 / 8.0 );       // inner edge dof
            WALBERLA_CHECK_FLOAT_EQUAL( it.value(), expected );
         }
      }
   }

   std::map< eo, uint_t > numNeighborElements = {
       { eo::X, 6 },
       { eo::Y, 4 },
       { eo::Z, 6 },
       { eo::XY, 6 },
       { eo::XZ, 4 },
       { eo::YZ, 6 },
       { eo::XYZ, 4 },

   };
   for ( auto it : FunctionIterator< EdgeDoFFunction< real_t > >( u.getEdgeDoFFunction(), lowerLevel ) )
   {
      if ( it.isOnMacroCell() )
      {
         if ( it.isEdgeDoF() )
         {
            real_t expected = 0;
            expected += 1.0;             // vertex
            expected += 2.0 * 3.0 / 4.0; // edgedofs on edge
            expected += real_c( numNeighborElements[it.edgeDoFOrientation()] ) * 2.0 *
                        ( 1.0 / 2.0 ); // edgedofs at edge (but different orientation)
            expected += real_c( numNeighborElements[it.edgeDoFOrientation()] ) * 1.0 *
                        ( 1.0 / 4.0 ); // edgedofs same orientation but not at edge
            expected += real_c( numNeighborElements[it.edgeDoFOrientation()] ) * ( 1.0 / 4.0 ); // inner edge dof
            WALBERLA_LOG_INFO_ON_ROOT( "Difference at edgedof (actual - expected): " << it.value() - expected )
            WALBERLA_LOG_INFO_ON_ROOT( it )
            WALBERLA_CHECK_FLOAT_EQUAL( it.value(), expected );
         }
      }
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   testWeightsInCell( 3 );
   testWeightsInCell( 4 );

   return 0;
}
