/*
 * Copyright (c) 2024 Nils Kohl.
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

#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

/// This test simply checks whether each DoF is traversed exactly once when interpolate() is called by counting the
/// number of invocations and comparing that number with the number of global DoFs.
template < typename FunctionType >
void runTest( const MeshInfo& meshInfo, const uint_t& level )
{
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, 1 );

   FunctionType u( "u", storage, level, level );

   uint_t numberOfInvocations = 0;

   std::function< real_t( const Point3D& x ) > countingFunction = [&]( const Point3D& x ) {
      numberOfInvocations++;
      return 1.0;
   };

   u.interpolate( countingFunction, level );

   uint_t numberOfDoFs = numberOfGlobalDoFs( u, level );

   WALBERLA_LOG_INFO_ON_ROOT( " ------- interpolate invocations: " << numberOfInvocations );
   WALBERLA_LOG_INFO_ON_ROOT( " ------- number of DoFs:          " << numberOfDoFs );

   WALBERLA_CHECK_EQUAL( numberOfInvocations, numberOfDoFs );

   bool vtk = true;

   if ( vtk )
   {
      auto vtkOutput = hyteg::VTKOutput( ".", "FunctionInterpolateOnceTest", storage );
      vtkOutput.add( u );
      vtkOutput.write( level );
   }
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::vector< MeshInfo > meshInfos;

   meshInfos.push_back( MeshInfo::fromGmshFile( "../../../data/meshes/tri_1el.msh" ) );
   meshInfos.push_back( MeshInfo::fromGmshFile( "../../../data/meshes/bfs_126el.msh" ) );
   meshInfos.push_back( MeshInfo::fromGmshFile( "../../../data/meshes/3D/tet_1el.msh" ) );
   meshInfos.push_back( MeshInfo::fromGmshFile( "../../../data/meshes/3D/cube_24el.msh" ) );

   uint_t minLevel = 0;
   uint_t maxLevel = 3;

   for ( const auto& meshInfo : meshInfos )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " - new mesh" )
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " --- level " << level )

         WALBERLA_LOG_INFO_ON_ROOT( " ----- P0 " )
         runTest< P0Function< real_t > >( meshInfo, level );

         WALBERLA_LOG_INFO_ON_ROOT( " ----- P1 " )
         runTest< P1Function< real_t > >( meshInfo, level );

         WALBERLA_LOG_INFO_ON_ROOT( " ----- P2 " )
         runTest< P2Function< real_t > >( meshInfo, level );
      }
   }

   return 0;
}