/*
* Copyright (c) 2017-2024 Marcus Mohr.
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

// The MicroMeshPositionTest checks that the following free-functions from
// the micromesh namespace work correctly:
//
// - micromesh::microVertexPosition()
// - micromesh::microEdgeCenterPosition() <-- overload version for EdgeDoFIndices
// - micromesh::microFaceCenterPosition()
//
// The testing is only done in 2D so far, but for all currently supported cases, i.e.
//
// a) no blending or micromesh
// b) blending
// c) use of degree 1 micromesh
// d) use of degree 2 micromesh

#include "core/Environment.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

// this class describes the coordinates of the single macro-face used for testing and
// provides a method to compute the position of micro-vertices on the computational domain
struct TestTriangle
{
   constexpr static const real_t a      = static_cast< real_t >( 1.0 );
   constexpr static const real_t deltaX = static_cast< real_t >( 2.0 );
   constexpr static const real_t deltaY = static_cast< real_t >( 0.8 );

   static Point3D vertexFromIndex( const indexing::Index& idx, uint_t level )
   {
      real_t  hx = deltaX / real_c( levelinfo::num_microedges_per_edge( level ) );
      real_t  hy = deltaY / real_c( levelinfo::num_microedges_per_edge( level ) );
      Point3D aux{ a + real_c( idx.x() + idx.y() ) * hx, a + real_c( idx.y() ) * hy, real_c( 0 ) };
      return aux;
   }
};

// utility function to pretty-print coordinates of a node
std::ostream& operator<<( std::ostream& out, const Point3D& point )
{
   out << "( " << std::scientific << std::setprecision( 2 ) << std::showpos << point( 0 ) << ", " << point( 1 ) << " )";
   return out;
}

// utility function to check identity of two nodes node by comparing their coordinates
void checkNodeIdentity( const Point3D& first, const Point3D& second )
{
   for ( int k = 0; k < 3; ++k )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( first[k], second[k] );
   }
}

std::array< Point3D, 3 > getMicroFaceVertexCoordinatesOnComputationalDomain( uint_t                 level,
                                                                             const indexing::Index& microFaceIndex,
                                                                             facedof::FaceType      faceType )
{
   std::array< Point3D, 3 > microVertices;
   const auto               microVertexIndices = facedof::macroface::getMicroVerticesFromMicroFace( microFaceIndex, faceType );

#ifdef MICRO_MESH_POSITION_TEST_DEBUG
   std::cout << ".................\n"
             << microVertexIndices[0] << '\n'
             << microVertexIndices[1] << '\n'
             << microVertexIndices[2] << '\n'
             << "faceType is " << facedof::FaceTypeToStr.at( faceType ) << '\n'
             << "microFaceIndex = " << microFaceIndex << '\n'
             << ".................\n"
             << std::endl;
#endif

   microVertices[0] = TestTriangle::vertexFromIndex( microVertexIndices[0], level );
   microVertices[1] = TestTriangle::vertexFromIndex( microVertexIndices[1], level );
   microVertices[2] = TestTriangle::vertexFromIndex( microVertexIndices[2], level );

   return microVertices;
}

Point3D mapPointToPhysicalDomain( const Point3D&                             point,
                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                  const PrimitiveID&                         faceID )
{
   Point3D mapped;

   const auto micromesh = storage->getMicroMesh();

   if ( !micromesh )
   {
      const auto& geometryMap = storage->getFace( faceID )->getGeometryMap();
      geometryMap->evalF( point, mapped );
   }
   else
   {
      WALBERLA_ABORT( "Do use different version of mapPointToPhysicalDomain for MicroMesh case." );
   }

   return mapped;
}

// Note: This function only works correctly if point is a micro-face vertex (micromesh degree 1 and 2)
// or an edge midpoint (micromesh degree 2 ) !
Point3D mapPointToPhysicalDomain( const Point3D&                                                   point,
                                  const std::shared_ptr< PrimitiveStorage >&                       storage,
                                  std::vector< std::function< real_t( const hyteg::Point3D& ) > >& map )
{
   Point3D mapped;

   const auto micromesh = storage->getMicroMesh();

   if ( micromesh != nullptr )
   {
      mapped( 0 ) = map[0]( point );
      mapped( 1 ) = map[1]( point );
      mapped( 2 ) = real_c( 0 );
   }
   else
   {
      WALBERLA_ABORT( "Do use different version of mapPointToPhysicalDomain for blending case." );
   }

   return mapped;
}

void run2DtestOnLevel( uint_t level, bool beVerbose = true )
{
   WALBERLA_LOG_INFO_ON_ROOT( "=====================================" );
   WALBERLA_LOG_INFO_ON_ROOT( " MicroMeshPositionTest: 2D, level " << level );
   WALBERLA_LOG_INFO_ON_ROOT( "=====================================" );

   const real_t a      = real_c( 1.0 );
   const real_t deltaX = real_c( 2.0 );
   const real_t deltaY = real_c( 0.8 );

   MeshInfo meshInfo = MeshInfo::singleTriangle( Point2D{ a, a }, Point2D{ a + deltaX, a }, Point2D{ a + deltaX, a + deltaY } );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // -------------------------------------------------------------------------------------------

   // for the different micromesh application scenarios (currently maximally four)
   const uint_t                                                    numScenarios = 4;
   std::array< std::shared_ptr< PrimitiveStorage >, numScenarios > stores;
   std::array< std::string, numScenarios >                         scenarioLabel;

   // Scenario #0: computational = physical domain
   stores[0]        = std::make_shared< PrimitiveStorage >( setupStorage );
   scenarioLabel[0] = "[w/o blending or micromesh]";

   // Scenario #1: with use of a blending map
   PolarCoordsMap::setMap( setupStorage );
   stores[1]        = std::make_shared< PrimitiveStorage >( setupStorage );
   scenarioLabel[1] = "[w/ blending]";

   // define a "blending" function to use for the micromesh
   std::vector< std::function< real_t( const hyteg::Point3D& ) > > polarCoords = {
       []( const Point3D& x ) { return x[0] * std::cos( x[1] ); }, []( const Point3D& x ) { return x[0] * std::sin( x[1] ); } };

   // Scenario #2: with a linear isoparametric mapping
   IdentityMap::setMap( setupStorage );
   stores[2]                   = std::make_shared< PrimitiveStorage >( setupStorage );
   const auto microMeshDegree1 = std::make_shared< micromesh::MicroMesh >( stores[2], level, level, 1, 2 );
   micromesh::interpolateAndCommunicate( *microMeshDegree1, polarCoords, level );
   stores[2]->setMicroMesh( microMeshDegree1 );
   scenarioLabel[2] = "[w/ micromesh (degree 1)]";

   // Scenario #3: with a quadratic isoparametric mapping
   stores[3]                   = std::make_shared< PrimitiveStorage >( setupStorage );
   const auto microMeshDegree2 = std::make_shared< micromesh::MicroMesh >( stores[3], level, level, 2, 2 );
   micromesh::interpolateAndCommunicate( *microMeshDegree2, polarCoords, level );
   stores[3]->setMicroMesh( microMeshDegree2 );
   scenarioLabel[3] = "[w/ micromesh (degree 2)]";

   // -------------------------------------------------------------------------------------------

   WALBERLA_LOG_INFO_ON_ROOT( "- testing vertex positions in 2D" );

   for ( uint_t k = 0; k < numScenarios; ++k )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  * testing setting #" << k << " " << scenarioLabel[k] );

      const auto& storage = stores[k];

      auto allFaceIDs = storage->getFaceIDs();
      WALBERLA_CHECK_EQUAL( allFaceIDs.size(), 1 );
      PrimitiveID faceID{ allFaceIDs[0] };

      for ( const auto& it : vertexdof::macroface::Iterator( level, 0 ) )
      {
         Point3D node     = micromesh::microVertexPosition( storage, faceID, level, it );
         Point3D unmapped = TestTriangle::vertexFromIndex( it, level );
         Point3D control;

         if ( k < 2 )
         {
            control = mapPointToPhysicalDomain( unmapped, storage, faceID );
         }
         else
         {
            control = mapPointToPhysicalDomain( unmapped, storage, polarCoords );
         }

         if ( beVerbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "    node = " << node << ", control = " << control );
         }

         checkNodeIdentity( node, control );
      }
   }

   // -------------------------------------------------------------------------------------------

   WALBERLA_LOG_INFO_ON_ROOT( "- testing edge midpoint positions in 2D" );

   for ( uint_t k = 0; k < numScenarios; ++k )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  * testing setting #" << k << " " << scenarioLabel[k] );

      const auto& storage = stores[k];

      for ( const auto& faceType : facedof::allFaceTypes )
      {
         for ( const auto& microFaceIndex : facedof::macroface::Iterator( level, faceType, 0 ) )
         {
            auto allFaceIDs = storage->getFaceIDs();
            WALBERLA_CHECK_EQUAL( allFaceIDs.size(), 1 );
            PrimitiveID faceID{ allFaceIDs[0] };

            // obtain coordinates of vertices of current micro-face
            std::array< Point3D, 3 > coords =
                getMicroFaceVertexCoordinatesOnComputationalDomain( level, microFaceIndex, faceType );

            // compute edge midpoints on computational domain
            std::array< Point3D, 3 > edgeCenters;
            edgeCenters[0] = coords[0] + ( coords[1] - coords[0] ) * real_c( 0.5 );
            edgeCenters[1] = coords[0] + ( coords[2] - coords[0] ) * real_c( 0.5 );
            edgeCenters[2] = coords[1] + ( coords[2] - coords[1] ) * real_c( 0.5 );

            if ( k == 0 || k == 1 )
            {
               edgeCenters[0] = mapPointToPhysicalDomain( edgeCenters[0], storage, faceID );
               edgeCenters[1] = mapPointToPhysicalDomain( edgeCenters[1], storage, faceID );
               edgeCenters[2] = mapPointToPhysicalDomain( edgeCenters[2], storage, faceID );
            }
            else if ( k == 2 )
            {
               std::array< Point3D, 3 > newCoords;
               newCoords[0] = mapPointToPhysicalDomain( coords[0], storage, polarCoords );
               newCoords[1] = mapPointToPhysicalDomain( coords[1], storage, polarCoords );
               newCoords[2] = mapPointToPhysicalDomain( coords[2], storage, polarCoords );

               edgeCenters[0] = newCoords[0] + ( newCoords[1] - newCoords[0] ) * real_c( 0.5 );
               edgeCenters[1] = newCoords[0] + ( newCoords[2] - newCoords[0] ) * real_c( 0.5 );
               edgeCenters[2] = newCoords[1] + ( newCoords[2] - newCoords[1] ) * real_c( 0.5 );
            }
            else if ( k == 3 )
            {
               edgeCenters[0] = mapPointToPhysicalDomain( edgeCenters[0], storage, polarCoords );
               edgeCenters[1] = mapPointToPhysicalDomain( edgeCenters[1], storage, polarCoords );
               edgeCenters[2] = mapPointToPhysicalDomain( edgeCenters[2], storage, polarCoords );
            }

            const auto microVertexIndices = facedof::macroface::getMicroVerticesFromMicroFace( microFaceIndex, faceType );

            std::array< indexing::Index, 3 > microEdgeIndex;
            microEdgeIndex[0] = edgedof::calcEdgeDoFIndex( microVertexIndices[0], microVertexIndices[1] );
            microEdgeIndex[1] = edgedof::calcEdgeDoFIndex( microVertexIndices[0], microVertexIndices[2] );
            microEdgeIndex[2] = edgedof::calcEdgeDoFIndex( microVertexIndices[1], microVertexIndices[2] );

            std::array< edgedof::EdgeDoFOrientation, 3 > microEdgeOrientation;
            microEdgeOrientation[0] = edgedof::calcEdgeDoFOrientation( microVertexIndices[0], microVertexIndices[1] );
            microEdgeOrientation[1] = edgedof::calcEdgeDoFOrientation( microVertexIndices[0], microVertexIndices[2] );
            microEdgeOrientation[2] = edgedof::calcEdgeDoFOrientation( microVertexIndices[1], microVertexIndices[2] );

            for ( uint_t j = 0; j < 3; ++j )
            {
               Point3D node =
                   micromesh::microEdgeCenterPosition( storage, faceID, level, microEdgeIndex[j], microEdgeOrientation[j] );
               Point3D control = edgeCenters[j];
               if ( beVerbose )
               {
                  WALBERLA_LOG_INFO_ON_ROOT( "    node = " << node << ", control = " << control );
               }
               checkNodeIdentity( node, control );
            }
         }
      }
   }

   // -------------------------------------------------------------------------------------------

   WALBERLA_LOG_INFO_ON_ROOT( "- testing face center position in 2D" );

   for ( uint_t k = 0; k < numScenarios; ++k )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "  * testing setting #" << k << " " << scenarioLabel[k] );

      const auto& storage = stores[k];

      auto allFaceIDs = storage->getFaceIDs();
      WALBERLA_CHECK_EQUAL( allFaceIDs.size(), 1 );
      PrimitiveID faceID{ allFaceIDs[0] };

      for ( const auto& faceType : facedof::allFaceTypes )
      {
         for ( const auto& microFaceIndex : facedof::macroface::Iterator( level, faceType, 0 ) )
         {
            Point3D center = micromesh::microFaceCenterPosition( storage, faceID, level, microFaceIndex, faceType );
            Point3D control;

            // obtain coordinates of vertices of current micro-face
            std::array< Point3D, 3 > coords =
                getMicroFaceVertexCoordinatesOnComputationalDomain( level, microFaceIndex, faceType );

            if ( k == 0 || k == 1 )
            {
               Point3D aux = ( coords[0] + coords[1] + coords[2] ) / real_c( 3 );
               control     = mapPointToPhysicalDomain( aux, storage, faceID );
            }
            else if ( k == 2 )
            {
               std::array< Point3D, 3 > newCoords;
               newCoords[0] = mapPointToPhysicalDomain( coords[0], storage, polarCoords );
               newCoords[1] = mapPointToPhysicalDomain( coords[1], storage, polarCoords );
               newCoords[2] = mapPointToPhysicalDomain( coords[2], storage, polarCoords );
               control      = ( newCoords[0] + newCoords[1] + newCoords[2] ) / real_c( 3 );
            }
            else if ( k == 3 )
            {
               const auto microVertexIndices = facedof::macroface::getMicroVerticesFromMicroFace( microFaceIndex, faceType );

               std::array< indexing::Index, 3 > microEdgeIndex;
               microEdgeIndex[0] = edgedof::calcEdgeDoFIndex( microVertexIndices[0], microVertexIndices[1] );
               microEdgeIndex[1] = edgedof::calcEdgeDoFIndex( microVertexIndices[0], microVertexIndices[2] );
               microEdgeIndex[2] = edgedof::calcEdgeDoFIndex( microVertexIndices[1], microVertexIndices[2] );

               std::array< edgedof::EdgeDoFOrientation, 3 > microEdgeOrientation;
               microEdgeOrientation[0] = edgedof::calcEdgeDoFOrientation( microVertexIndices[0], microVertexIndices[1] );
               microEdgeOrientation[1] = edgedof::calcEdgeDoFOrientation( microVertexIndices[0], microVertexIndices[2] );
               microEdgeOrientation[2] = edgedof::calcEdgeDoFOrientation( microVertexIndices[1], microVertexIndices[2] );

               std::array< Point3D, 6 > newCoords;
               newCoords[0] = mapPointToPhysicalDomain( coords[0], storage, polarCoords );
               newCoords[1] = mapPointToPhysicalDomain( coords[1], storage, polarCoords );
               newCoords[2] = mapPointToPhysicalDomain( coords[2], storage, polarCoords );
               newCoords[3] =
                   micromesh::microEdgeCenterPosition( storage, faceID, level, microEdgeIndex[0], microEdgeOrientation[0] );
               newCoords[4] =
                   micromesh::microEdgeCenterPosition( storage, faceID, level, microEdgeIndex[1], microEdgeOrientation[1] );
               newCoords[5] =
                   micromesh::microEdgeCenterPosition( storage, faceID, level, microEdgeIndex[2], microEdgeOrientation[2] );
               control = -( newCoords[0] + newCoords[1] + newCoords[2] ) / real_c( 9 ) +
                         ( newCoords[3] + newCoords[4] + newCoords[5] ) * real_c( 4 ) / real_c( 9 );
            }

            if ( beVerbose )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "    center = " << center << ", control = " << control );
            }

            checkNodeIdentity( center, control );
         }
      }
   }
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   for ( uint_t level = 0; level <= 3; ++level )
   {
      run2DtestOnLevel( level, false );
   }

   return 0;
}
