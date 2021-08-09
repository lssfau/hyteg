/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Marcus Mohr.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

// Perform some basic test to check that methods of P[12]VectorFunctions
// can be instantiated, called, executed and exported

namespace hyteg {

template < typename vfType >
static void testVectorFunction( bool beVerbose, std::string tag, std::string typeName )
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   vfType vec_f( "vecFunc", storage, minLevel, maxLevel );
   vfType aux_f( "auxFunc", storage, minLevel, maxLevel );

   // Testing traits
   FunctionTrait< vfType > vfKind;
   WALBERLA_LOG_INFO_ON_ROOT( "Detected typename is: '" << vfKind.getTypeName() << "'" );
   WALBERLA_CHECK_EQUAL( typeName, FunctionTrait< vfType >::getTypeName() );

   // Interpolate
   std::function< real_t( const hyteg::Point3D& ) > xComp = []( const Point3D& ) { return real_c( 2 ); };
   std::function< real_t( const hyteg::Point3D& ) > yComp = []( const Point3D& x ) { return x[0] * x[1]; };

   walberla::WcTimingPool timer;

   timer["Interpolate"].start();
   vec_f.interpolate( { xComp, yComp }, maxLevel, DoFType::All );
   timer["Interpolate"].end();

   // Assign
   timer["Assign"].start();
   aux_f.assign( { 3.0 }, { vec_f }, maxLevel, DoFType::All );
   timer["Assign"].end();

   // Add
   timer["Add"].start();
   aux_f.add( { { 4.0, 3.0 } }, { { vec_f, vec_f } }, maxLevel, DoFType::All );
   timer["Add"].end();

   // GetMaxMagnitude
   timer["GetMaxComponentMagnitude"].start();
   aux_f.getMaxComponentMagnitude( maxLevel, DoFType::All );
   timer["GetMaxComponentMagnitude"].end();

   // Dot
   timer["Dot"].start();
   const real_t scalarProduct = aux_f.dotGlobal( vec_f, maxLevel, DoFType::All );
   timer["Dot"].end();
   WALBERLA_LOG_INFO_ON_ROOT( "dot product = " << scalarProduct );

   // Output VTK
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = tag + "VectorFunctionExportViaComponents";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( dynamic_cast< typename vfType::VectorComponentType& >( vec_f[0] ) );
      vtkOutput.add( dynamic_cast< typename vfType::VectorComponentType& >( vec_f[1] ) );
      vtkOutput.write( maxLevel );

      std::string fName2 = tag + "VectorFunctionExport";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName2 << "'" );
      VTKOutput vtkOutput2( fPath, fName2, storage );
      vtkOutput2.add( vec_f );
      vtkOutput2.write( maxLevel );
   }

   WALBERLA_LOG_INFO_ON_ROOT( timer );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing P1VectorFunction" );
   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );
   hyteg::testVectorFunction< hyteg::P1VectorFunction< walberla::real_t > >( true, "P1", "P1VectorFunction" );

   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );
   WALBERLA_LOG_INFO_ON_ROOT( " Testing P2VectorFunction" );
   WALBERLA_LOG_INFO_ON_ROOT( "==========================" );
   hyteg::testVectorFunction< hyteg::P2VectorFunction< walberla::real_t > >( true, "P2", "P2VectorFunction" );

   return EXIT_SUCCESS;
}
