/*
 * Copyright (c) 2017-2025 Marcus Mohr.
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
#include "hyteg/experimental/P2PlusBubbleFunction.hpp"

#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

void testInterpolation( std::shared_ptr< PrimitiveStorage > storage, std::array< real_t, 3 > tolerance, bool doVTKOutput = false )
{
   const uint_t minLevel = 4;
   const uint_t maxLevel = minLevel;

   auto function0 = P2PlusBubbleFunction< real_t >( "poly0", storage, minLevel, maxLevel );
   auto function1 = P2PlusBubbleFunction< real_t >( "poly1", storage, minLevel, maxLevel );
   auto function2 = P2PlusBubbleFunction< real_t >( "poly2", storage, minLevel, maxLevel );
   auto function3 = P2PlusBubbleFunction< real_t >( "poly3", storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > polyDegree0 = []( const Point3D& ) -> real_t { return real_c( 1 ); };

   std::function< real_t( const hyteg::Point3D& ) > polyDegree1 = []( const Point3D& x ) -> real_t {
      return real_c( 2 ) * x[1] - x[0];
   };

   std::function< real_t( const hyteg::Point3D& ) > polyDegree2 = []( const Point3D& x ) -> real_t {
      return x[0] * x[0] + real_c( 2 ) * x[0] * x[1];
   };

   std::function< real_t( const hyteg::Point3D& ) > polyDegree3 = []( const Point3D& x ) -> real_t {
      return x[0] * x[1] * ( real_c( 1 ) - x[0] - x[1] ) * real_c( 27 );
   };

   function0.interpolate( polyDegree0, maxLevel, All );
   function1.interpolate( polyDegree1, maxLevel, All );
   function2.interpolate( polyDegree2, maxLevel, All );
   function3.interpolate( polyDegree3, maxLevel, All );

   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "p2_plus_bubble", storage );
      vtkOutput.add( function0 );
      vtkOutput.add( function1 );
      vtkOutput.add( function2 );
      vtkOutput.add( function3 );
      vtkOutput.write( maxLevel );
   }

   // for the quadratic polynomial all bubble dofs need to be zero
   const auto& bubbleFunc0 = function0.getVolumeDoFFunction();
   const auto& bubbleFunc1 = function1.getVolumeDoFFunction();
   const auto& bubbleFunc2 = function2.getVolumeDoFFunction();
   const auto& bubbleFunc3 = function3.getVolumeDoFFunction();

   real_t checkValue0 = bubbleFunc0.getMaxDoFMagnitude( maxLevel );
   real_t checkValue1 = bubbleFunc1.getMaxDoFMagnitude( maxLevel );
   real_t checkValue2 = bubbleFunc2.getMaxDoFMagnitude( maxLevel );
   real_t checkValue3 = bubbleFunc3.getMaxDoFMagnitude( maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 0: max magnitude of bubble dofs = " << checkValue0 )
   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 1: max magnitude of bubble dofs = " << checkValue1 )
   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 2: max magnitude of bubble dofs = " << checkValue2 )
   WALBERLA_LOG_INFO_ON_ROOT( " * polynomial 3: max magnitude of bubble dofs = " << checkValue3 )

   WALBERLA_CHECK_FLOAT_EQUAL( checkValue0, real_c( 0 ) );
   WALBERLA_CHECK_LESS( checkValue1, tolerance[0] );
   WALBERLA_CHECK_LESS( checkValue2, tolerance[1] );
   WALBERLA_CHECK_LESS( checkValue3, tolerance[2] );
}

void runInterpolationTests()
{
   // simple mesh
   MeshInfo              mesh1 = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/tri_1el.msh" ) );
   SetupPrimitiveStorage setupStorage1( mesh1, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // annulus
   MeshInfo              mesh2 = MeshInfo::meshAnnulus( real_c( 1 ), real_c( 2 ), MeshInfo::CROSS, 4, 2 );
   SetupPrimitiveStorage setupStorage2( mesh2, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorage2 );

   // unit square with affine blending
   MeshInfo mesh3 = MeshInfo::meshRectangle(
       Point2D( real_c( -1 ), real_c( -1 ) ), Point2D( real_c( +1 ), real_c( +1 ) ), MeshInfo::CRISSCROSS, 1, 1 );
   SetupPrimitiveStorage setupStorage3( mesh3, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // define our affine map (rotation + scaling + shift)
   Matrix2r mat;
   real_t   phi = real_c( 2.0 / 9.0 ) * walberla::math::pi;
   mat( 0, 0 )  = +std::cos( phi );
   mat( 0, 1 )  = -std::sin( phi );
   mat( 1, 0 )  = +std::sin( phi ) * real_c( 2.25 );
   mat( 1, 1 )  = +std::cos( phi ) * real_c( 2.25 );
   Point2D vec( real_c( -7.0 ), real_c( 3.0 ) );
   AffineMap2D::setMap( setupStorage3, mat, vec );

   // quadratic mapping
   MeshInfo mesh4 = MeshInfo::meshRectangle(
       Point2D( real_c( -1 ), real_c( -1 ) ), Point2D( real_c( +1 ), real_c( +1 ) ), MeshInfo::CRISSCROSS, 1, 1 );
   SetupPrimitiveStorage setupStorage4( mesh3, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::vector< std::function< real_t( const hyteg::Point3D& ) > > mapper = {
       []( const Point3D& x ) { return x[0]; },
       []( const Point3D& x ) {
          return x[1] + ( real_c( 1 ) - x[0] * x[0] ) * x[1] * ( std::sqrt( real_c( 2 ) ) - real_c( 1 ) );
       } };
   std::shared_ptr< PrimitiveStorage > storage4 = std::make_shared< PrimitiveStorage >( setupStorage4, 1 );
   uint_t                              level    = 4;
   const auto microMeshDegree2                  = std::make_shared< micromesh::MicroMesh >( storage4, level, level, 2, 2 );
   micromesh::interpolateAndCommunicate( *microMeshDegree2, mapper, level );
   storage4->setMicroMesh( microMeshDegree2 );

   std::shared_ptr< PrimitiveStorage > storage1 = std::make_shared< PrimitiveStorage >( setupStorage1, 1 );
   std::shared_ptr< PrimitiveStorage > storage2 = std::make_shared< PrimitiveStorage >( setupStorage2, 1 );
   std::shared_ptr< PrimitiveStorage > storage3 = std::make_shared< PrimitiveStorage >( setupStorage3, 1 );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation on unmodified mesh" );
   testInterpolation( storage1, { real_c( 3e-16 ), real_c( 3e-16 ), real_c( 2.5e-4 ) } );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation on mesh with affine blending" );
   testInterpolation( storage3, { real_c( 8e-15 ), real_c( 3e-14 ), real_c( 2.1e-3 ) } );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation with quadratic micromesh" );
   testInterpolation( storage4, { real_c( 9e-16 ), real_c( 4.5e-5 ), real_c( 3e-3 ) } );

   WALBERLA_LOG_INFO_ON_ROOT( "Testing interpolation on mesh with annulus blending" );
   testInterpolation( storage2, { real_c( 1.8e-4 ), real_c( 5e-4 ), real_c( 0.03 ) } );
}

void runDiffusionTest( bool doVTKOutput = false )
{
   // setup domain as unit square (0,1)^2
   MeshInfo mesh = MeshInfo::meshRectangle(
       Point2D( real_c( 0 ), real_c( 0 ) ), Point2D( real_c( +1 ), real_c( +1 ) ), MeshInfo::CRISSCROSS, 1, 1 );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // pick an analytical solution for the homogeneous diffusion equation
   real_t freqNum{ real_c( 2 ) };

   std::function< real_t( const hyteg::Point3D& ) > expr_analytic = [&freqNum]( const Point3D& x ) {
      using walberla::math::pi;
      return std::sinh( freqNum * pi * x[1] ) / std::sinh( freqNum * pi ) * std::sin( freqNum * pi * x[0] );
   };

   // setup some P2+ functions
   uint_t minLevel = 0;
   uint_t maxLevel = 4;

   P2PlusBubbleFunction< real_t > u_analytic( "true solution", storage, minLevel, maxLevel );
   u_analytic.interpolate( expr_analytic, maxLevel, All );

   P2PlusBubbleFunction< real_t > weak_diffusion( "weak diffusion", storage, minLevel, maxLevel );
   // weak_diffusion.interpolate( expr_analytic, maxLevel, All );

   // output stuff for visualisation
   if ( doVTKOutput )
   {
      VTKOutput vtkOutput( ".", "p2+_diffusion", storage );
      vtkOutput.add( u_analytic );
      vtkOutput.write( maxLevel );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   using namespace hyteg;

   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   runInterpolationTests();

   runDiffusionTest( true );

   return EXIT_SUCCESS;
}
