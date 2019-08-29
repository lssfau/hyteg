
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/FunctionIterator.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

void testWeightsInCell( const uint_t& lowerLevel )
{
   typedef edgedof::EdgeDoFOrientation eo;

   const auto storage = PrimitiveStorage::createFromGmshFile( "../../data/meshes/3D/tet_1el.msh" );

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
       {eo::X, 6},
       {eo::Y, 4},
       {eo::Z, 6},
       {eo::XY, 6},
       {eo::XZ, 4},
       {eo::YZ, 6},
       {eo::XYZ, 4},

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