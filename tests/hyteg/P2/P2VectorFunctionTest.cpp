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

#include "hyteg/p2functionspace/P2VectorFunction.hpp"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

// Perform some basic test to check that methods of P2VectorFunction
// can be called and executed (and instantiated, see restrict and
// prolongate below ;-)

namespace hyteg {

static void testP2Function()
{
   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   MeshInfo                            mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P2VectorFunction< real_t > vec_f( "vecFunc", storage, minLevel, maxLevel );
   P2VectorFunction< real_t > aux_f( "auxFunc", storage, minLevel, maxLevel );

   // Interpolate
   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& ) -> real_t { return real_c( 2 ); };

   walberla::WcTimingPool timer;

   timer["Interpolate"].start();
   vec_f.interpolate( expr, maxLevel, DoFType::All );
   timer["Interpolate"].end();

   // Assign
   timer["Assign"].start();
   aux_f.assign( {3.0}, {vec_f}, maxLevel, DoFType::All );
   timer["Assign"].end();

   // Add
   timer["Add"].start();
   aux_f.add( {{4.0, 3.0}}, {{vec_f, vec_f}}, maxLevel, DoFType::All );
   timer["Add"].end();

   // Dot
   timer["Dot"].start();
   const real_t scalarProduct = aux_f.dotGlobal( vec_f, maxLevel, DoFType::All );
   timer["Dot"].end();
   WALBERLA_LOG_INFO_ON_ROOT( "dot product = " << scalarProduct );

   // Output VTK
   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "P2VectorFunctionTest";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( vec_f[0] );
      vtkOutput.add( vec_f[1] );
      vtkOutput.write( maxLevel );

      std::string fName2 = "P2VectorFunctionExport";
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
   hyteg::testP2Function();

   return EXIT_SUCCESS;
}
