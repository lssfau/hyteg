/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr.
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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

// Perform basic I/O test with VTKOutput objects

// using namespace hyteg::n1e1;

namespace hyteg {

static void exportFunctions2D( uint_t level )
{
   uint_t minLevel = level;
   uint_t maxLevel = level;

   // MeshInfo                         mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/penta_5el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage   = std::make_shared< PrimitiveStorage >( setupStorage );
   std::shared_ptr< PrimitiveStorage > storageDG = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // Setup some functions
   P0Function< real_t > p0ScalarFunc1( "P0 scalar function 1", storageDG, minLevel, maxLevel );
   P0Function< real_t > p0ScalarFunc2( "P0 scalar function 2", storageDG, minLevel, maxLevel );

   P1Function< real_t > p1ScalarFunc1( "P1 scalar function 1", storage, minLevel, maxLevel );
   P1Function< real_t > p1ScalarFunc2( "P1 scalar function 2", storage, minLevel, maxLevel );

   P2Function< real_t > p2ScalarFunc1( "P2 scalar function 1", storage, minLevel, maxLevel );
   P2Function< real_t > p2ScalarFunc2( "P2 scalar function 2", storage, minLevel, maxLevel );

   DG1Function< real_t > dg1ScalarFunc1( "DG1 scalar function 1", storageDG, minLevel, maxLevel );
   DG1Function< real_t > dg1ScalarFunc2( "DG1 scalar function 2", storageDG, minLevel, maxLevel );

   P1VectorFunction< real_t > p1VectorFunc( "P1 vector function", storage, minLevel, maxLevel );
   P2VectorFunction< real_t > p2VectorFunc( "P2 vector function", storage, minLevel, maxLevel );
   EGFunction< real_t >       egVectorFunc( "EG vector function", storageDG, minLevel, maxLevel );

   // Interpolate
   p0ScalarFunc1.interpolate( 1.0, maxLevel, DoFType::All );
   p0ScalarFunc2.interpolate( 2.0, maxLevel, DoFType::All );

   p1ScalarFunc1.interpolate( 1.0, maxLevel, DoFType::All );
   p1ScalarFunc2.interpolate( 2.0, maxLevel, DoFType::All );

   std::function< real_t( const hyteg::Point3D& ) > xFunc = []( const Point3D& p ) -> real_t { return -2.0 * p[0]; };
   std::function< real_t( const hyteg::Point3D& ) > yFunc = []( const Point3D& p ) -> real_t { return p[0] + p[1]; };
   std::vector< std::function< real_t( const hyteg::Point3D& ) > > vecExpr = { xFunc, yFunc };
   p1VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );
   p2VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );
   egVectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );

   // Output VTK
   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "VTKOutputTest-P0";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput0( fPath, fName, storageDG );
      vtkOutput0.add( p0ScalarFunc1 );
      vtkOutput0.add( p0ScalarFunc2 );
      vtkOutput0.write( maxLevel );

      fName = "VTKOutputTest-P1";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput1( fPath, fName, storage );
      vtkOutput1.add( p1ScalarFunc1 );
      vtkOutput1.add( p1ScalarFunc2 );
      vtkOutput1.add( p1VectorFunc );
      vtkOutput1.write( maxLevel );

      fName = "VTKOutputTest-P2";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput2( fPath, fName, storage );
      vtkOutput2.add( p2ScalarFunc1 );
      vtkOutput2.add( p2ScalarFunc2 );
      vtkOutput2.add( p2VectorFunc );
      vtkOutput2.write( maxLevel );

      fName = "VTKOutputTest-DG1";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutputDG1( fPath, fName, storageDG );
      vtkOutputDG1.add( dg1ScalarFunc1 );
      vtkOutputDG1.add( dg1ScalarFunc2 );
      vtkOutputDG1.write( maxLevel );

      // Do EG functions and also check "CRTP inheritance"
      fName = "VTKOutputTest-EG";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      FEFunctionWriter< VTKOutput >* vtkOutputEG = new VTKOutput( fPath, fName, storageDG );
      vtkOutputEG->add( egVectorFunc );
      vtkOutputEG->write( maxLevel );
   }
}

static void exportFunctions3D( uint_t level )
{
   uint_t minLevel = level;
   uint_t maxLevel = level;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage   = std::make_shared< PrimitiveStorage >( setupStorage );
   std::shared_ptr< PrimitiveStorage > storageDG = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // Expressions for interpolation
   // std::function< real_t( const hyteg::Point3D& ) > xFunc = []( const Point3D& p ) -> real_t { return -2.0*p[0]; };
   // std::function< real_t( const hyteg::Point3D& ) > yFunc = []( const Point3D& p ) -> real_t { return p[0]+p[1]+p[2]; };
   // std::function< real_t( const hyteg::Point3D& ) > zFunc = []( const Point3D& p ) -> real_t { return 3.0*p[0]+p[2]; };
   std::function< real_t( const hyteg::Point3D& ) >                xFunc   = []( const Point3D& p ) -> real_t { return p[0]; };
   std::function< real_t( const hyteg::Point3D& ) >                yFunc   = []( const Point3D& p ) -> real_t { return p[1]; };
   std::function< real_t( const hyteg::Point3D& ) >                zFunc   = []( const Point3D& p ) -> real_t { return p[2]; };
   std::vector< std::function< real_t( const hyteg::Point3D& ) > > vecExpr = { xFunc, yFunc, zFunc };
   std::function< Point3D( const Point3D& ) >                      vecFunc = []( const Point3D& p ) {
      Point3D result{ p[0], p[1], p[2] };
      return result;
   };

   // Setup some functions
   P0Function< real_t > p0ScalarFunc1( "P0 scalar function 1", storageDG, minLevel, maxLevel );
   P0Function< real_t > p0ScalarFunc2( "P0 scalar function 2", storageDG, minLevel, maxLevel );
   P0Function< real_t > p0ScalarFunc3( "P0 scalar function 3", storageDG, minLevel, maxLevel );

   P1Function< real_t > p1ScalarFunc1( "P1 scalar function 1", storage, minLevel, maxLevel );
   P1Function< real_t > p1ScalarFunc2( "P1 scalar function 2", storage, minLevel, maxLevel );
   P1Function< real_t > p1ScalarFunc3( "P1 scalar function 3", storage, minLevel, maxLevel );

   P2Function< real_t > p2ScalarFunc1( "P2 scalar function 1", storage, minLevel, maxLevel );
   P2Function< real_t > p2ScalarFunc2( "P2 scalar function 2", storage, minLevel, maxLevel );
   P1Function< real_t > p2ScalarFunc3( "P2 scalar function 3", storage, minLevel, maxLevel );

   DG1Function< real_t > dg1ScalarFunc1( "DG1 scalar function 1", storageDG, minLevel, maxLevel );
   DG1Function< real_t > dg1ScalarFunc2( "DG1 scalar function 2", storageDG, minLevel, maxLevel );
   DG1Function< real_t > dg1ScalarFunc3( "DG1 scalar function 3", storageDG, minLevel, maxLevel );

   P1VectorFunction< real_t > p1VectorFunc( "P1 vector function", storage, minLevel, maxLevel );
   P2VectorFunction< real_t > p2VectorFunc( "P2 vector function", storage, minLevel, maxLevel );

   EGFunction< real_t > egVectorFunc( "EG vector function", storageDG, minLevel, maxLevel );

   n1e1::N1E1VectorFunction< real_t > n1e1VectorFunc( "N1E1 vector function", storage, minLevel, maxLevel );

   // Interpolate
   p0ScalarFunc1.interpolate( vecExpr[0], maxLevel, DoFType::All );
   p0ScalarFunc2.interpolate( vecExpr[1], maxLevel, DoFType::All );
   p0ScalarFunc3.interpolate( vecExpr[2], maxLevel, DoFType::All );

   p1ScalarFunc1.interpolate( vecExpr[0], maxLevel, DoFType::All );
   p1ScalarFunc2.interpolate( vecExpr[1], maxLevel, DoFType::All );
   p1ScalarFunc3.interpolate( vecExpr[2], maxLevel, DoFType::All );

   p2ScalarFunc1.interpolate( vecExpr[0], maxLevel, DoFType::All );
   p2ScalarFunc2.interpolate( vecExpr[1], maxLevel, DoFType::All );
   p2ScalarFunc3.interpolate( vecExpr[2], maxLevel, DoFType::All );

   p1VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );
   p2VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );
   egVectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );

   n1e1VectorFunc.interpolate( vecFunc, maxLevel, DoFType::All );

   // Output VTK
   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "VTKOutputTest3D-P0";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutputP0( fPath, fName, storageDG );
      vtkOutputP0.add( p0ScalarFunc1 );
      vtkOutputP0.add( p0ScalarFunc2 );
      vtkOutputP0.add( p0ScalarFunc3 );
      vtkOutputP0.setVTKDataFormat( vtk::DataFormat::ASCII );
      vtkOutputP0.write( maxLevel );

      fName = "VTKOutputTest3D-P1";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutputP1( fPath, fName, storage );
      vtkOutputP1.add( p1ScalarFunc1 );
      vtkOutputP1.add( p1ScalarFunc2 );
      vtkOutputP1.add( p1ScalarFunc3 );
      vtkOutputP1.add( p1VectorFunc );
      vtkOutputP1.setParameter( "vtkDataFormat", "BINARY" );
      vtkOutputP1.write( maxLevel );

      fName = "VTKOutputTest3D-P2";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutputP2( fPath, fName, storage );
      vtkOutputP2.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutputP2.add( p2ScalarFunc1 );
      vtkOutputP2.add( p2ScalarFunc2 );
      vtkOutputP2.add( p2ScalarFunc3 );
      vtkOutputP2.add( p2VectorFunc );
      vtkOutputP2.write( maxLevel );

      fName = "VTKOutputTest3D-DG1";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutputDG1( fPath, fName, storageDG );
      vtkOutputDG1.add( dg1ScalarFunc1 );
      vtkOutputDG1.add( dg1ScalarFunc2 );
      vtkOutputDG1.add( dg1ScalarFunc3 );
      vtkOutputDG1.write( maxLevel );

      fName = "VTKOutputTest3D-N1E1";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutputN1E1( fPath, fName, storage );
      vtkOutputN1E1.add( n1e1VectorFunc );
      vtkOutputN1E1.write( maxLevel );

      fName = "VTKOutputTest3D-EG";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutputEG( fPath, fName, storageDG );
      vtkOutputEG.add( egVectorFunc );
      // would currently fail as 3D export of EGFunction is not fully implemented, yet
      vtkOutputEG.write( maxLevel );
   }
}

template < typename value_t >
static void exportIntegerFunctions()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 2;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/penta_5el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // Setup some functions
   P0Function< value_t >       p0Enumerator( "P0", storage, minLevel, maxLevel );
   P1Function< value_t >       p1Enumerator( "P1", storage, minLevel, maxLevel );
   P2Function< value_t >       p2Enumerator( "P2", storage, minLevel, maxLevel );
   EdgeDoFFunction< value_t >  edEnumerator( "EdgeDoF", storage, minLevel, maxLevel );
   P2VectorFunction< value_t > v2Enumerator( "P2Vector", storage, minLevel, maxLevel );

   // Fill with values
   p0Enumerator.enumerate( maxLevel );
   p1Enumerator.enumerate( maxLevel );
   p2Enumerator.enumerate( maxLevel );
   edEnumerator.enumerate( maxLevel );
   v2Enumerator.enumerate( maxLevel );

   // Output VTK
   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "VTKOutputTestIntFuncs" + walberla::vtk::typeToString< value_t >();
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( p0Enumerator );
      vtkOutput.add( p1Enumerator );
      vtkOutput.add( p2Enumerator );
      vtkOutput.add( edEnumerator );
      vtkOutput.write( maxLevel );
   }
}

static void testVTKQuadraticTriangle( uint_t meshType, uint_t level )
{
   const uint_t minLevel = level;
   const uint_t maxLevel = level;

   using walberla::math::pi;

   std::shared_ptr< MeshInfo > mesh;

   switch ( meshType )
   {
   case 1:
      mesh = std::make_shared< MeshInfo >( MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ) );
      break;
   case 2:
      mesh = std::make_shared< MeshInfo >( MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ) );
      break;
   case 3:
      mesh = std::make_shared< MeshInfo >( MeshInfo::fromGmshFile( "../../data/meshes/penta_5el.msh" ) );
      break;
   case 4:
      mesh = std::make_shared< MeshInfo >( MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ) );
      break;
   case 5:
      mesh = std::make_shared< MeshInfo >(
          MeshInfo::meshAnnulus( real_c( 1 ), real_c( 2 ), real_c( 0.25 ) * pi, real_c( 0.75 ) * pi, MeshInfo::CROSS, 4, 2 ) );
      break;
   default:
      WALBERLA_ABORT( "meshType = " << meshType << " is not supported!" );
   }

   SetupPrimitiveStorage setupStorage( *mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   if ( meshType == 3 )
      AnnulusMap::setMap( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Setup test function
   P2Function< real_t >                             p2ScalarFunc( "P2 scalar function", storage, minLevel, maxLevel );
   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& p ) -> real_t {
      return real_c( -2.0 ) * p[0] + p[1] + real_c( 3 );
   };
   p2ScalarFunc.interpolate( expr, maxLevel, All );

   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "VTKQuadraticTriangle";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( p2ScalarFunc );
      vtkOutput.write( maxLevel );
   }
}

static void testVTKQuadraticTetra( uint_t meshType, uint_t level )
{
   uint_t minLevel = level;
   uint_t maxLevel = level;

   using walberla::math::pi;

   std::shared_ptr< MeshInfo > mesh;

   switch ( meshType )
   {
   case 1:
      mesh = std::make_shared< MeshInfo >( MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ) );
      break;
   case 2:
      mesh = std::make_shared< MeshInfo >( MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ) );
      break;
   case 3:
      mesh = std::make_shared< MeshInfo >( MeshInfo::fromGmshFile( "../../data/meshes/3D/three_tets_with_two_joint_faces.msh" ) );
      break;
   default:
      WALBERLA_ABORT( "meshType = " << meshType << " is not supported!" );
   }

   SetupPrimitiveStorage               setupStorage( *mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Setup test function
   P2Function< real_t >                             p2ScalarFunc( "P2 scalar function", storage, minLevel, maxLevel );
   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& p ) -> real_t {
      return real_c( -2.0 ) * p[0] + p[1] + real_c( 3 ) + p[2] * p[0];
   };
   p2ScalarFunc.interpolate( expr, maxLevel, All );

   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "VTKQuadraticTetra";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutput.add( p2ScalarFunc );
      vtkOutput.write( maxLevel );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for 2D meshes:" );
   for ( uint_t level = 0; level <= 2; ++level )
   {
      hyteg::exportFunctions2D( level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for 3D meshes:" );
   for ( uint_t level = 0; level <= 2; ++level )
   {
      hyteg::exportFunctions3D( level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for value_t = int32_t:" );
   hyteg::exportIntegerFunctions< int32_t >();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for value_t = int64_t:" );
   hyteg::exportIntegerFunctions< int64_t >();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export with VTK_QUADRATIC_TRIANGLE:" );
   for ( uint_t level = 0; level <= 2; ++level )
   {
      hyteg::testVTKQuadraticTriangle( 2, level );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export with VTK_QUADRATIC_TETRA:" );
   for ( uint_t level = 0; level <= 2; ++level )
   {
      hyteg::testVTKQuadraticTetra( 3, level );
   }

   return EXIT_SUCCESS;
}
