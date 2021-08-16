/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Marcus Mohr.
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

#include "hyteg/dataexport/VTKOutput.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

// Perform basic I/O test with VTKOutput objects

namespace hyteg {

static void exportFunctions2D()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 2;

   // MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/penta_5el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Setup some functions
   P1Function< real_t > p1ScalarFunc1( "P1 scalar function 1", storage, minLevel, maxLevel );
   P1Function< real_t > p1ScalarFunc2( "P1 scalar function 2", storage, minLevel, maxLevel );

   P2Function< real_t > p2ScalarFunc1( "P2 scalar function 1", storage, minLevel, maxLevel );
   P2Function< real_t > p2ScalarFunc2( "P2 scalar function 2", storage, minLevel, maxLevel );

   P1VectorFunction< real_t > p1VectorFunc( "P1 vector function", storage, minLevel, maxLevel );
   P2VectorFunction< real_t > p2VectorFunc( "P2 vector function", storage, minLevel, maxLevel );

   // Interpolate
   p1ScalarFunc1.interpolate( 1.0, maxLevel, DoFType::All );
   p1ScalarFunc2.interpolate( 2.0, maxLevel, DoFType::All );

   std::function< real_t( const hyteg::Point3D& ) > xFunc = []( const Point3D& p ) -> real_t { return -2.0 * p[0]; };
   std::function< real_t( const hyteg::Point3D& ) > yFunc = []( const Point3D& p ) -> real_t { return p[0] + p[1]; };
   std::vector< std::function< real_t( const hyteg::Point3D& ) > > vecExpr = { xFunc, yFunc };
   p1VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );
   p2VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );

   // Output VTK
   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "VTKOutputTest";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( p1ScalarFunc1 );
      vtkOutput.add( p1ScalarFunc2 );
      vtkOutput.add( p1VectorFunc );
      vtkOutput.write( maxLevel );

      fName = "VTKOutputTestP2";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput2( fPath, fName, storage );
      vtkOutput2.add( p2ScalarFunc1 );
      vtkOutput2.add( p2ScalarFunc2 );
      vtkOutput2.add( p2VectorFunc );
      vtkOutput2.write( maxLevel );
   }
}

static void exportFunctions3D()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 2;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Expressions for interpolation
   // std::function< real_t( const hyteg::Point3D& ) > xFunc = []( const Point3D& p ) -> real_t { return -2.0*p[0]; };
   // std::function< real_t( const hyteg::Point3D& ) > yFunc = []( const Point3D& p ) -> real_t { return p[0]+p[1]+p[2]; };
   // std::function< real_t( const hyteg::Point3D& ) > zFunc = []( const Point3D& p ) -> real_t { return 3.0*p[0]+p[2]; };
   std::function< real_t( const hyteg::Point3D& ) >                xFunc   = []( const Point3D& p ) -> real_t { return p[0]; };
   std::function< real_t( const hyteg::Point3D& ) >                yFunc   = []( const Point3D& p ) -> real_t { return p[1]; };
   std::function< real_t( const hyteg::Point3D& ) >                zFunc   = []( const Point3D& p ) -> real_t { return p[2]; };
   std::vector< std::function< real_t( const hyteg::Point3D& ) > > vecExpr = { xFunc, yFunc, zFunc };

   // Setup some functions
   P1Function< real_t > p1ScalarFunc1( "P1 scalar function 1", storage, minLevel, maxLevel );
   P1Function< real_t > p1ScalarFunc2( "P1 scalar function 2", storage, minLevel, maxLevel );
   P1Function< real_t > p1ScalarFunc3( "P1 scalar function 3", storage, minLevel, maxLevel );

   P2Function< real_t > p2ScalarFunc1( "P2 scalar function 1", storage, minLevel, maxLevel );
   P2Function< real_t > p2ScalarFunc2( "P2 scalar function 2", storage, minLevel, maxLevel );
   P1Function< real_t > p2ScalarFunc3( "P2 scalar function 3", storage, minLevel, maxLevel );

   P1VectorFunction< real_t > p1VectorFunc( "P1 vector function", storage, minLevel, maxLevel );
   P2VectorFunction< real_t > p2VectorFunc( "P2 vector function", storage, minLevel, maxLevel );

   // Interpolate
   p1ScalarFunc1.interpolate( vecExpr[0], maxLevel, DoFType::All );
   p1ScalarFunc2.interpolate( vecExpr[1], maxLevel, DoFType::All );
   p1ScalarFunc3.interpolate( vecExpr[2], maxLevel, DoFType::All );

   p2ScalarFunc1.interpolate( vecExpr[0], maxLevel, DoFType::All );
   p2ScalarFunc2.interpolate( vecExpr[1], maxLevel, DoFType::All );
   p2ScalarFunc3.interpolate( vecExpr[2], maxLevel, DoFType::All );

   p1VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );
   p2VectorFunc.interpolate( vecExpr, maxLevel, DoFType::All );

   // Output VTK
   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "VTKOutputTest3D";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( p1ScalarFunc1 );
      vtkOutput.add( p1ScalarFunc2 );
      vtkOutput.add( p1ScalarFunc3 );
      vtkOutput.add( p1VectorFunc );
      vtkOutput.write( maxLevel );

      fName = "VTKOutputTest3DP2";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput2( fPath, fName, storage );
      vtkOutput2.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutput2.add( p2ScalarFunc1 );
      vtkOutput2.add( p2ScalarFunc2 );
      vtkOutput2.add( p2ScalarFunc3 );
      vtkOutput2.add( p2VectorFunc );
      vtkOutput2.write( maxLevel );
   }
}

template < typename value_t >
static void exportIntegerFunctions()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 2;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/penta_5el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Setup some functions
   P1Function< value_t >       p1Enumerator( "P1", storage, minLevel, maxLevel );
   P2Function< value_t >       p2Enumerator( "P2", storage, minLevel, maxLevel );
   EdgeDoFFunction< value_t >  edEnumerator( "EdgeDoF", storage, minLevel, maxLevel );
   DGFunction< value_t >       dgEnumerator( "DGDoF", storage, minLevel, maxLevel );
   P2VectorFunction< value_t > v2Enumerator( "P2Vector", storage, minLevel, maxLevel );

   // Fill with values
   p1Enumerator.enumerate( maxLevel );
   p2Enumerator.enumerate( maxLevel );
   dgEnumerator.enumerate( maxLevel );
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
      vtkOutput.add( p1Enumerator );
      vtkOutput.add( p2Enumerator );
      vtkOutput.add( dgEnumerator );
      vtkOutput.add( edEnumerator );
      vtkOutput.add( v2Enumerator );
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
   hyteg::exportFunctions2D();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for 3D meshes:" );
   hyteg::exportFunctions3D();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for value_t = int32_t:" );
   hyteg::exportIntegerFunctions< int32_t >();

   WALBERLA_LOG_INFO_ON_ROOT( "Testing export for value_t = int64_t:" );
   hyteg::exportIntegerFunctions< int64_t >();

   return EXIT_SUCCESS;
}
