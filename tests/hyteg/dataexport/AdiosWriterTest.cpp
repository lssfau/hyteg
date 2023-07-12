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
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"

#include <adios2.h>
#include <cstring>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// This test checks that exporting FE functions with the AdiosWriter
// class compiles and runs through without aborting.

using namespace hyteg;

template < typename func_t >
void initScalarFunction( const func_t& function, const uint_t level )
{
   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& p ) -> real_t {
      return real_c( -2 ) * p[0] + real_c( 3 ) * p[1];
   };
   function.interpolate( expr, level, DoFType::All );
}

template < typename func_t >
void initVectorFunction( const func_t& function, const uint_t level )
{
   std::function< real_t( const hyteg::Point3D& ) > exprX = []( const Point3D& p ) -> real_t {
      return real_c( -2 ) * p[0] + real_c( 3 ) * p[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > exprY = []( const Point3D& p ) -> real_t {
      return std::sin( p[0] ) * std::cos( real_c( 3 ) * p[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > exprZ = []( const Point3D& p ) -> real_t {
      return std::sin( p[0] + p[2] ) * std::cos( real_c( 3 ) * p[1] * p[2] );
   };
   function.interpolate( { exprX, exprY, exprZ }, level, DoFType::All );
}

template < typename value_t >
void runTest( std::shared_ptr< PrimitiveStorage > storage, std::string baseFileName, uint_t level )
{
   P1Function< real_t > p1Func( "P1TestFunction", storage, level, level );
   initScalarFunction( p1Func, level );
   P1VectorFunction< real_t > p1VecFunc( "P1VectorTestFunction", storage, level, level );
   initVectorFunction( p1VecFunc, level );

   P2Function< real_t > p2Func( "P2TestFunction", storage, level, level );
   initScalarFunction( p2Func, level );
   P2VectorFunction< real_t > p2VecFunc( "P2VectorTestFunction", storage, level, level );
   initVectorFunction( p2VecFunc, level );

   P2P1TaylorHoodFunction< real_t > stokesFunc( "StokesTestFunction", storage, level, level );

   AdiosWriter adiosWriter( "../../output", baseFileName, storage );
   adiosWriter.add( p1Func );
   adiosWriter.add( p2Func );
   adiosWriter.add( p1VecFunc );
   adiosWriter.add( p2VecFunc );
   adiosWriter.add( stokesFunc );

   // perform some write steps
   adiosWriter.write( level, 0 );

   // should produce a warning
   P1Function< real_t > foo( "foo", storage, level, level );
   adiosWriter.add( foo );
   
   adiosWriter.write( level, 1 );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Execute test with a 2D mesh
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing AdiosWriter for 2D Mesh ***" );
   {
      MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/LShape_6el.msh" );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
      runTest< real_t >( storage, "AdiosWriterTest_2D", 3 );
   }

   // Execute test with a 3D mesh
   WALBERLA_LOG_INFO_ON_ROOT( "*** Testing AdiosWriter for 3D Mesh ***" );
   {
      MeshInfo              mesh = MeshInfo::meshSphericalShell( 2, 2, real_c( 1 ), real_c( 2 ) );
      SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      IcosahedralShellMap::setMap( setupStorage );
      setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
      runTest< real_t >( storage, "AdiosWriterTest_3D", 2 );
   }

   return EXIT_SUCCESS;
}
