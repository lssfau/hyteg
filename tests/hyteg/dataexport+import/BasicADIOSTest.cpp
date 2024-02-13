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
#include <adios2.h>
#include <cstring>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

// This test sets up a P1Function on a single triangle macro mesh and exports the
// complete data held by the macro face via ADIOS2. It then reads the data in again
// to setup a cloned function and compares it to the original. This is done twice
// once for the BP5 engine (so file IO) and then with the "inline" engine.

using namespace hyteg;

void initSrcFunction( const P1Function< real_t >& srcFunc, const uint_t level )
{
   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& p ) -> real_t {
      return real_c( -2 ) * p[0] + real_c( 3 ) * p[1];
   };
   srcFunc.interpolate( expr, level, DoFType::All );
}

void compareFunctions( const P1Function< real_t >& original, const P1Function< real_t >& clone, uint_t level, bool writeOutVTK )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Comparing cloned function to orignal:" );

   auto                 storage = original.getStorage();
   P1Function< real_t > check( "auxilliary function", storage, level, level );
   check.assign( { real_c( 1 ), real_c( -1 ) }, { original, clone }, level, All );
   real_t measure = check.getMaxMagnitude( level, All );

   if ( writeOutVTK )
   {
      VTKOutput vtkOutput( ".", "check", storage );
      vtkOutput.add( original );
      vtkOutput.add( clone );
      vtkOutput.add( check );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_FLOAT_EQUAL( measure, real_c( 0 ) );

   WALBERLA_LOG_INFO_ON_ROOT( "-> checks out" );
}

// test cloning a function via file export / import
void testFileDataExchange( const P1Function< real_t >&          original,
                           std::shared_ptr< PrimitiveStorage >& storage,
                           const uint_t                         level,
                           adios2::ADIOS&                       adios,
                           bool                                 writeOutVTK = false )
{
   std::string filename{ "BasicAdiosTest.bp" };

   P1Function< real_t > clone( "cloned function", storage, level, level );

   adios2::IO bpWriter = adios.DeclareIO( "BPWriter" );
   adios2::IO bpReader = adios.DeclareIO( "BPReader" );

   bpWriter.SetEngine( "BP5" );
   bpReader.SetEngine( "BP5" );

   uint_t nDoFs = numberOfLocalDoFs< P1FunctionTag >( *( original.getStorage() ), level );
   WALBERLA_LOG_INFO( "Going to write " << nDoFs << " real_t values into '" << filename << "'" );

   {
      // make sure face halo is up-to-date
      communication::syncFunctionBetweenPrimitives( original, level );

      // define a variable to represent the actual data
      adios2::Variable< real_t > varFaceData = bpWriter.DefineVariable< real_t >( "Face Data", {}, {}, { nDoFs } );

      // define an attribute to store the function name
#if 0
      adios2::Attribute< std::string > funcName = bpWriter.DefineAttribute( "Function Name", original.getFunctionName(), "Face Data" );
      // attribute is associated to the variable, so we do not explicitely export it;
      // which means the C++ variable is unused :(
      WALBERLA_UNUSED( funcName );
#else
      // or we can discard the return value
      bpWriter.DefineAttribute( "Function Name", original.getFunctionName(), "Face Data" );
#endif

      auto    face     = storage->getFaces().begin()->second;
      real_t* faceData = face->getData( original.getFaceDataID() )->getPointer( level );

      adios2::Engine engineOutput = bpWriter.Open( filename, adios2::Mode::Write );
      engineOutput.BeginStep();
      engineOutput.Put( varFaceData, faceData );
      engineOutput.EndStep();
      engineOutput.Close();
   }

   {
      auto    face     = storage->getFaces().begin()->second;
      real_t* faceData = face->getData( clone.getFaceDataID() )->getPointer( level );

      adios2::Engine engineInput = bpReader.Open( filename, adios2::Mode::Read );
      engineInput.BeginStep();

      adios2::Variable< real_t > varFaceData = bpReader.InquireVariable< real_t >( "Face Data" );
      if ( !varFaceData )
      {
         WALBERLA_ABORT( "Could not find 'Face Data' in '" << filename << "'" );
      }
      else
      {
         WALBERLA_LOG_INFO( "Going to read values back in from '" << filename << "'" );
      }

      adios2::Attribute< std::string > funcName = bpReader.InquireAttribute< std::string >( "Function Name", "Face Data" );
      WALBERLA_LOG_INFO_ON_ROOT( "Attribute name is " << funcName.Name() );
      WALBERLA_LOG_INFO_ON_ROOT( "Attribute type is " << funcName.Type() );
      std::string origName = funcName.Data()[0];
      WALBERLA_LOG_INFO_ON_ROOT( "Cloning function '" << origName << "'" );
      WALBERLA_CHECK( origName == original.getFunctionName() );

      engineInput.Get( varFaceData, faceData );
      engineInput.EndStep();
      engineInput.Close();

      // propagate info from face halo to lower-dim primitives
      clone.communicateAdditively< Face, Edge >( level, true );
      clone.communicateAdditively< Face, Vertex >( level, true );
   }

   // check that clone matches original
   compareFunctions( original, clone, level, writeOutVTK );
}

// test cloning via "inline" data exchange
void testInlineDataExchange( const P1Function< real_t >&          original,
                             std::shared_ptr< PrimitiveStorage >& storage,
                             const uint_t                         level,
                             adios2::ADIOS&                       adios,
                             bool                                 writeOutVTK = false )
{
   P1Function< real_t > clone( "cloned function", storage, level, level );
   P1Function< real_t > check( "auxilliary function", storage, level, level );

   // reader and writer share an IO object
   adios2::IO inlineIO = adios.DeclareIO( "InlineIO" );
   inlineIO.SetEngine( "Inline" );

   uint_t                     nDoFs       = numberOfLocalDoFs< P1FunctionTag >( *( original.getStorage() ), level );
   adios2::Variable< real_t > varFaceData = inlineIO.DefineVariable< real_t >( "Face Data", {}, {}, { nDoFs } );

   // make sure face halo is up-to-date
   communication::syncFunctionBetweenPrimitives( original, level );

   // need to open both reader and writer before starting IO in inline mode
   adios2::Engine inlineWriter = inlineIO.Open( "inlinePut", adios2::Mode::Write );
   adios2::Engine inlineReader = inlineIO.Open( "inlineFetch", adios2::Mode::Read );

   auto    face       = storage->getFaces().begin()->second;
   real_t* dataOrig   = face->getData( original.getFaceDataID() )->getPointer( level );
   real_t* dataClone  = face->getData( clone.getFaceDataID() )->getPointer( level );
   real_t* dataBuffer = nullptr;

   WALBERLA_LOG_INFO( "Going to put " << nDoFs << " real_t values" );
   inlineWriter.BeginStep();
   inlineWriter.Put( varFaceData, dataOrig );
   inlineWriter.EndStep();

   WALBERLA_LOG_INFO( "Going to fetch values" );
   inlineReader.BeginStep();
   inlineReader.Get( varFaceData, &dataBuffer );
   // it is, of course, a little awkward to now do a copy ;-) but I see no easy way
   // to replace the actual data buffer of the face
   std::memcpy( dataClone, dataBuffer, sizeof( real_t ) * nDoFs );
   inlineReader.EndStep();

   // close engines again
   inlineReader.Close();
   inlineWriter.Close();

   // propagate info from face halo to lower-dim primitives
   clone.communicateAdditively< Face, Edge >( level, true );
   clone.communicateAdditively< Face, Vertex >( level, true );

   // check that clone matches original
   compareFunctions( original, clone, level, writeOutVTK );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo              mesh = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t level = 5;

   P1Function< real_t > original( "original function", storage, level, level );
   initSrcFunction( original, level );

#ifdef WALBERLA_BUILD_WITH_MPI
   adios2::ADIOS adios( walberla::MPIManager::instance()->comm() );
#else
   adios2::ADIOS adios;
#endif

   WALBERLA_LOG_INFO_ON_ROOT( "*** Test Cloning via FileIO ***" );
   testFileDataExchange( original, storage, level, adios );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "*** Test Inline Cloning ***" );
   testInlineDataExchange( original, storage, level, adios, true );

   return EXIT_SUCCESS;
}
