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
// test that the product one^T*M*one with mass matrix M and vector of ones gives area of domain
#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// Hunting the NaN
#include <cfenv>

using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

template < typename OperatorType >
void checkArea( std::shared_ptr< PrimitiveStorage > storage,
                real_t                              area,
                std::string                         tag,
                const uint_t                        minLevel  = 2,
                real_t                              tolerance = -1.0,
                bool                                outputVTK = false )
{
   // const uint_t minLevel = 2;
   const uint_t maxLevel = storage->hasGlobalCells() ? 3 : 4;

   OperatorType massOp( storage, minLevel, maxLevel );

   typename OperatorType::srcType aux( "aux", storage, minLevel, maxLevel );
   typename OperatorType::srcType vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );

   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      vecOfOnes.interpolate( real_c( 1.0 ), lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dotGlobal( aux, lvl );

      if ( outputVTK )
      {
         VTKOutput vtkOutput( "../../output", tag, storage );
         vtkOutput.add( vecOfOnes );
         vtkOutput.add( aux );
         vtkOutput.write( lvl );
      }

      if ( tolerance < 0.0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure << " (" << tag << ")" );
         WALBERLA_CHECK_FLOAT_EQUAL( measure, area );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure << ", difference = " << std::abs( measure - area )
                                                 << " (" << tag << ")" );
         WALBERLA_CHECK_LESS( std::abs( measure - area ), tolerance );
      }
   }
}

void setMap( SetupPrimitiveStorage& setupStorage, std::shared_ptr< GeometryMap > map )
{
   for ( auto it : setupStorage.getFaces() )
   {
      setupStorage.setGeometryMap( it.second->getID(), map );
   }
   for ( auto it : setupStorage.getEdges() )
   {
      setupStorage.setGeometryMap( it.second->getID(), map );
   }
   for ( auto it : setupStorage.getVertices() )
   {
      setupStorage.setGeometryMap( it.second->getID(), map );
   }
}

void logSectionHeader( const char* header )
{
   std::string hdr( header );
   size_t      len = hdr.length();
   std::string separator( len + 2, '-' );
   WALBERLA_LOG_INFO_ON_ROOT( separator << "\n " << hdr << "\n" << separator );
}

int main( int argc, char** argv )
{
#ifndef __APPLE__ 
    #ifndef _MSC_VER
        feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
   #endif
#endif
   walberla::debug::enterTestMode();

   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::unique_ptr< SetupPrimitiveStorage > setStore;
   std::shared_ptr< PrimitiveStorage >      primStore;

   Matrix2r mat;
   Point2D vec;

   // ----------
   //  2D Tests
   // ----------

   // Test with rectangle
   logSectionHeader( "Testing with RECTANGLE" );
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0.0, -1.0} ), Point2D( {2.0, 3.0} ), MeshInfo::CRISSCROSS, 1, 2 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1ConstantMassOperator >( primStore, 8.0, "P1ConstantMassOperator" );
   checkArea< P2ConstantMassOperator >( primStore, 8.0, "P2ConstantMassOperator" );
   checkArea< P1ElementwiseMassOperator >( primStore, 8.0, "P1ElementwiseMassOperator" );
   checkArea< P2ElementwiseMassOperator >( primStore, 8.0, "P2ElementwiseMassOperator" );

   // Test with backward facing step
   logSectionHeader( "Testing with BFS" );
   meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1ConstantMassOperator >( primStore, 1.75, "P1ConstantMassOperator" );
   checkArea< P2ConstantMassOperator >( primStore, 1.75, "P2ConstantMassOperator" );
   checkArea< P1ElementwiseMassOperator >( primStore, 1.75, "P1ElementwiseMassOperator" );
   checkArea< P2ElementwiseMassOperator >( primStore, 1.75, "P2ElementwiseMassOperator" );

   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 1.75, "P2ElementwiseBlendingMassOperator" );

   // ----------
   //  3D Tests
   // ----------

   // test with cuboid
   logSectionHeader( "Testing with Cuboid" );
   meshInfo = MeshInfo::meshCuboid( Point3D( {-1.0, -1.0, 0.0} ), Point3D( {2.0, 0.0, 2.0} ), 1, 2, 1 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1ConstantMassOperator >( primStore, 6.0, "P1ConstantMassOperator" );
   checkArea< P2ConstantMassOperator >( primStore, 6.0, "P2ConstantMassOperator" );
   checkArea< P1ElementwiseMassOperator >( primStore, 6.0, "P1ElementwiseMassOperator" );
   checkArea< P2ElementwiseMassOperator >( primStore, 6.0, "P2ElementwiseMassOperator" );

   // Test with coarse representation of thick spherical shell
   logSectionHeader( "Testing with Icosahedral Shell" );
   meshInfo = MeshInfo::meshSphericalShell( 2, {1.0, 2.0} );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore         = std::make_shared< PrimitiveStorage >( *setStore.get() );
   real_t edgeLength = 8.0 / ( std::sqrt( 10.0 + 2.0 * std::sqrt( 5.0 ) ) );
   real_t volume     = 5.0 / 12.0 * ( 3.0 + std::sqrt( 5.0 ) ) * edgeLength * edgeLength * edgeLength;
   edgeLength /= 2.0;
   volume -= 5.0 / 12.0 * ( 3.0 + std::sqrt( 5.0 ) ) * edgeLength * edgeLength * edgeLength;

   checkArea< P1ConstantMassOperator >( primStore, volume, "P1ConstantMassOperator" );
   checkArea< P2ConstantMassOperator >( primStore, volume, "P2ConstantMassOperator" );
   checkArea< P1ElementwiseMassOperator >( primStore, volume, "P1ElementwiseMassOperator" );
   checkArea< P2ElementwiseMassOperator >( primStore, volume, "P2ElementwiseMassOperator" );

   // -------------------
   //  2D Blending Tests
   // -------------------

   // Test with annulus
   logSectionHeader( "Testing with BLENDING( ANNULUS -- PolarCoordsMap )" );
   meshInfo = MeshInfo::meshRectangle( Point2D( {1.0, 0.0} ), Point2D( {2.0, 2.0 * pi} ), MeshInfo::CROSS, 1, 6 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setMap( *setStore.get(), std::make_shared< PolarCoordsMap >() );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1BlendingMassOperator >( primStore, 3.0 * pi, "P1BlendingMassOperator" );
   checkArea< P1ElementwiseBlendingMassOperator >( primStore, 3.0 * pi, "P1ElementwiseBlendingMassOperator" );
   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 3.0 * pi, "P2ElementwiseBlendingMassOperator" );

   // Test with annulus v2 (why do we need more refinement, compared to PolarCoordsMap?)
   logSectionHeader( "Testing with BLENDING( ANNULUS -- AnnulusMap )" );
   meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, 0.0, 2.0 * pi, MeshInfo::CROSS, 24, 4 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( *setStore.get() );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1ElementwiseBlendingMassOperator >( primStore, 3.0 * pi, "P1ElementwiseBlendingMassOperator", 3 );
   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 3.0 * pi, "P2ElementwiseBlendingMassOperator", 3 );
   checkArea< P1BlendingMassOperator >( primStore, 3.0 * pi, "P1BlendingMassOperator", 3 );

   // Test with unit square containing circular hole
   logSectionHeader( "Testing with BLENDING( SQUARE with CIRCULAR HOLE )" );
   meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/unitsquare_with_circular_hole.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{{0.5, 0.5, 0}};
   real_t  circleRadius = 0.25;

   for ( const auto& it : setStore.get()->getFaces() )
   {
      Face& face = *( it.second );

      std::vector< PrimitiveID > neighborEdgesOnBoundary = face.neighborEdges();
      neighborEdgesOnBoundary.erase(
          std::remove_if( neighborEdgesOnBoundary.begin(),
                          neighborEdgesOnBoundary.end(),
                          [&setStore]( const PrimitiveID& id ) { return !setStore.get()->onBoundary( id ); } ),
          neighborEdgesOnBoundary.end() );

      if ( neighborEdgesOnBoundary.size() > 0 )
      {
         Edge& edge = *( setStore.get()->getEdge( neighborEdgesOnBoundary[0] ) );

         if ( ( edge.getCoordinates()[0] - circleCenter ).norm() < 0.4 )
         {
            setStore.get()->setGeometryMap(
                edge.getID(), std::make_shared< CircularMap >( face, *setStore.get(), circleCenter, circleRadius ) );
            setStore.get()->setGeometryMap(
                face.getID(), std::make_shared< CircularMap >( face, *setStore.get(), circleCenter, circleRadius ) );
         }
      }
   }
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1BlendingMassOperator >( primStore, 1.0 - pi / 16.0, "P1BlendingMassOperator", 3 );
   checkArea< P1ElementwiseBlendingMassOperator >( primStore, 1.0 - pi / 16.0, "P1ElementwiseBlendingMassOperator", 3 );
   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 1.0 - pi / 16.0, "P2ElementwiseBlendingMassOperator", 3 );

   // Test with backward facing step and affine mapping
   logSectionHeader( "Testing with BLENDING( AFFINE_MAP rotation )" );
   meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   vec = Point2D( {0.0, 0.0} );
   real_t phi = real_c(27) * pi / real_c(180);
   mat(0,0) = +std::cos(phi);
   mat(0,1) = -std::sin(phi);
   mat(1,0) = +std::sin(phi);
   mat(1,1) = +std::cos(phi);
   AffineMap2D::setMap( *setStore.get(), mat, vec );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1BlendingMassOperator >( primStore, 1.75, "P1BlendingMassOperator" );
   checkArea< P1ElementwiseBlendingMassOperator >( primStore, 1.75, "P1ElementwiseBlendingMassOperator" );
   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 1.75, "P2ElementwiseBlendingMassOperator" );

   // Test with backward facing step and affine mapping
   logSectionHeader( "Testing with BLENDING( AFFINE_MAP shear, scale + shift )" );
   meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   vec = Point2D( {-2.0, 3.0} );
   real_t scalFac = real_c(2);
   mat(0,0) = real_c(scalFac);
   mat(0,1) = real_c(1);
   mat(1,0) = real_c(0);
   mat(1,1) = real_c(scalFac);
   AffineMap2D::setMap( *setStore.get(), mat, vec );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P1BlendingMassOperator >( primStore, 1.75*scalFac*scalFac, "P1BlendingMassOperator", 2, -1.0, true );
   checkArea< P1ElementwiseBlendingMassOperator >( primStore, 1.75*scalFac*scalFac, "P1ElementwiseBlendingMassOperator" );
   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 1.75*scalFac*scalFac, "P2ElementwiseBlendingMassOperator" );

   // -------------------
   //  3D Blending Tests
   // -------------------

   // Test with identity mapping
   logSectionHeader( "Testing with BLENDING( UNIT CUBE with IdentityMap )" );

   Point3D lowerLeftFront( {0.0, 0.0, 0.0} );
   Point3D upperRightBack( {1.0, 1.0, 1.0} );
   meshInfo = MeshInfo::meshCuboid( lowerLeftFront, upperRightBack, 1, 1, 1 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 1.0, "P2ElementwiseBlendingMassOperator", 2, 6e-8 );

   // Test with affine mapping
   logSectionHeader( "Testing with BLENDING( UNIT CUBE with AffineMap )" );

   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // define our affine map
   Matrix3r matAffineMap;
   matAffineMap(0,0) = +8.660254037844387e-01;
   matAffineMap(0,1) = -1.545084971874737e-01;
   matAffineMap(0,2) = -9.510565162951534e-01;
   matAffineMap(1,0) = +0.000000000000000e+00;
   matAffineMap(1,1) = +9.510565162951535e-01;
   matAffineMap(1,2) = -6.180339887498948e-01;
   matAffineMap(2,0) = +4.999999999999999e-01;
   matAffineMap(2,1) = +2.676165673298174e-01;
   matAffineMap(2,2) = +1.647278207092664e+00;
   Point3D vecAffineMap( {-7.0, 3.0, 2.0} );
   AffineMap3D::setMap( *setStore.get(), matAffineMap, vecAffineMap );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P2ElementwiseBlendingMassOperator >( primStore, 2.0, "P2ElementwiseBlendingMassOperator", 2, 6e-8 );

   // Test with thick spherical shell
   logSectionHeader( "Testing with BLENDING( Thick Spherical Shell -- IcosahedralShellMap )" );
   meshInfo = MeshInfo::meshSphericalShell( 2, 2, 1.0, 2.0 );
   setStore =
       std::make_unique< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   IcosahedralShellMap::setMap( *setStore.get() );
   primStore = std::make_shared< PrimitiveStorage >( *setStore.get() );

   checkArea< P2ElementwiseBlendingMassOperator >(
       primStore, 4.0 / 3.0 * pi * 7.0, "P2ElementwiseBlendingMassOperator", 2, 5e-6 );

   return EXIT_SUCCESS;
}
