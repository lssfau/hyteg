/*
 * Copyright (c) 2025 Ponsuganth Ilangovan
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

#include <core/Environment.h>
#include <core/math/Constants.h>

#include "hyteg/MeshQuality.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo meshInfo = MeshInfo::meshSphericalShell(3u, 2u, 1.22, 2.22);
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );

   const std::string storageFilename{ "sph_storage_mesh_file" };

   WALBERLA_ROOT_SECTION()
   {
      setupStorage->writeToFile( storageFilename, 1u );
#ifdef HYTEG_BUILD_WITH_ADIOS2
      setupStorage->writeToFile( walberla::format("%s.bp", storageFilename.c_str()), 1u, true );
#endif
   }

   // Load storage from file
   auto storageRead = std::make_shared< hyteg::PrimitiveStorage >( storageFilename, 1u );


   auto info     = storage->getGlobalInfo();
   auto infoRead = storageRead->getGlobalInfo();

   WALBERLA_LOG_INFO_ON_ROOT( "Storage built during run time:" )
   WALBERLA_LOG_INFO_ON_ROOT( info );

   WALBERLA_LOG_INFO_ON_ROOT( "Storage read from file:" )
   WALBERLA_LOG_INFO_ON_ROOT( infoRead );


   WALBERLA_CHECK_EQUAL( info, infoRead );

   const uint_t level = 1u;

   P2Function< real_t > TWithStorage("TWithStorage", storage, level, level);
   P2Function< real_t > TWithStorageRead("TWithStorageRead", storageRead, level, level);

   P2VectorFunction< real_t > uWithStorage("TWithStorage", storage, level, level);
   P2VectorFunction< real_t > uWithStorageRead("TWithStorageRead", storageRead, level, level);

   MMOCTransport< P2Function< real_t > > mmocTransportWithStorage(storage, level, level, TimeSteppingScheme::RK4);
   MMOCTransport< P2Function< real_t > > mmocTransportWithStorageRead(storageRead, level, level, TimeSteppingScheme::RK4);
   
   std::function< real_t(const Point3D&) > someFunc = [](const Point3D& x)
   {
      return std::sin(x[0]) * x[1] * x[1] * x[0] + std::cos(x[1]) * x[0] * x[0];
   };

   TWithStorage.interpolate(someFunc, level, All);
   TWithStorageRead.interpolate(someFunc, level, All);

   uWithStorage.interpolate(someFunc, level, All);
   uWithStorageRead.interpolate(someFunc, level, All);

   real_t uMax = uWithStorage.getMaxComponentMagnitude(level, All);
   real_t hMin = MeshQuality::getMinimalEdgeLength(storage, level);

   real_t dt = uMax / hMin;

   WALBERLA_LOG_INFO_ON_ROOT("dt = " << dt);

   mmocTransportWithStorage.step(TWithStorage, uWithStorage, uWithStorage, level, All, dt, 1u);
   mmocTransportWithStorageRead.step(TWithStorageRead, uWithStorageRead, uWithStorageRead, level, All, dt, 1u);

   TWithStorageRead.assign({1.0, -1.0}, {TWithStorageRead, TWithStorage}, level, All);

   real_t TWithStorageReadNorm = TWithStorageRead.dotGlobal(TWithStorageRead, level, All);

   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TWithStorageReadNorm       = %4.15e", TWithStorageReadNorm ) );
   // WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TWithStorageReadAdios2Norm = %4.15e", TWithStorageReadAdios2Norm ) );

   WALBERLA_CHECK_EQUAL(TWithStorageReadNorm, 0.0);

#ifdef HYTEG_BUILD_WITH_ADIOS2   
   // Load storage from adios2 file
   auto storageReadAdios2 = std::make_shared< hyteg::PrimitiveStorage >( walberla::format("%s.bp", storageFilename.c_str()), 1u, true );
   auto infoReadAdios2 = storageReadAdios2->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( "Storage read from Adios2 file:" )
   WALBERLA_LOG_INFO_ON_ROOT( infoReadAdios2 );
   WALBERLA_CHECK_EQUAL( info, infoReadAdios2 );

   P2Function< real_t > TWithStorageReadAdios2("TWithStorageReadAdios2", storageReadAdios2, level, level);
   P2VectorFunction< real_t > uWithStorageReadAdios2("TWithStorageReadAdios2", storageReadAdios2, level, level);
   MMOCTransport< P2Function< real_t > > mmocTransportWithStorageReadAdios2(storageReadAdios2, level, level, TimeSteppingScheme::RK4);

   TWithStorageReadAdios2.interpolate(someFunc, level, All);
   uWithStorageReadAdios2.interpolate(someFunc, level, All);

   mmocTransportWithStorageReadAdios2.step(TWithStorageReadAdios2, uWithStorageReadAdios2, uWithStorageReadAdios2, level, All, dt, 1u);
   
   TWithStorageReadAdios2.assign({1.0, -1.0}, {TWithStorageReadAdios2, TWithStorage}, level, All);

   real_t TWithStorageReadAdios2Norm = TWithStorageReadAdios2.dotGlobal(TWithStorageReadAdios2, level, All);
   WALBERLA_CHECK_EQUAL(TWithStorageReadAdios2Norm, 0.0);
#endif

   return EXIT_SUCCESS;
}
