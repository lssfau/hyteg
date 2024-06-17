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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dataexport/VTKOutput/VTKPointCloudOutput.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"

#include "terraneo/helpers/RadialProfiles.hpp"

using namespace hyteg;

namespace terraneo {

std::shared_ptr< PrimitiveStorage >
    setupSphericalShellStorage( const uint_t& nTan, const uint_t& nRad, const real_t& rMin, const real_t& rMax )
{
   auto meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax ) );

   SetupPrimitiveStorage setupStorage( *meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   IcosahedralShellMap::setMap( setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   return storage;
}

std::shared_ptr< PrimitiveStorage >
    setupSphericalShellStorage( const uint_t& nTan, std::vector<real_t> layers )
{
   auto meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshSphericalShell( nTan, layers ) );

   SetupPrimitiveStorage setupStorage( *meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   IcosahedralShellMap::setMap( setupStorage );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   return storage;
}

void testShellMath()
{
   // For the thick spherical shell, nRad is the number of radial SHELLS, not layers.
   // So the number of shells on refinement level 0 must be the same as nRad.
   // Considering P1 functions here (vertices).
   for ( uint_t nRad = 1; nRad < 10; nRad++ )
   {
      WALBERLA_CHECK_EQUAL( numberOfShells( nRad, 0, 1 ), nRad );
   }
}

template < typename FunctionType >
void testRadialPointCloudOutput( uint_t nTan, uint_t nRad, real_t rMin, real_t rMax, uint_t level )
{
   auto storage = setupSphericalShellStorage( nTan, nRad, rMin, rMax );

   FunctionType u( "u", storage, level, level );

   auto radius = []( const Point3D& x ) { return x.norm(); };

   u.interpolate( radius, level );

   RadialShellData< FunctionType > shellData;

   shellData.addDataFromFunction( u, rMin, rMax, nRad, level );

   uint_t numPointsOnAllShellsCombined = 0;
   for ( uint_t shell = 0; shell < numberOfShells( nRad, level, polynomialDegreeOfBasisFunctions< FunctionType >() ); shell++ )
   {
      auto numPointsOnShell = shellData.points( shell ).size();
      walberla::mpi::allReduceInplace( numPointsOnShell, walberla::mpi::Operation::SUM );
      WALBERLA_LOG_INFO_ON_ROOT( "shell " << shell << " | points: " << numPointsOnShell );

      numPointsOnAllShellsCombined += numPointsOnShell;
   }
   // That number should be the same as the number of DoFs!
   auto numGlobalDoFs = numberOfGlobalDoFs( u, level );

   const bool testableFunctionType = std::is_same_v< typename FunctionType::Tag, P1FunctionTag > ||
                                     std::is_same_v< typename FunctionType::Tag, P2FunctionTag > ||
                                     std::is_same_v< typename FunctionType::Tag, P1VectorFunctionTag > ||
                                     std::is_same_v< typename FunctionType::Tag, P2VectorFunctionTag >;
   WALBERLA_CHECK( testableFunctionType, "The next WALBERLA_CHECK_EQUAL() is not applicable for all function types :)" )
   WALBERLA_CHECK_EQUAL( numPointsOnAllShellsCombined, numGlobalDoFs / u.getDimension() );

   // Write one file per shell.
   for ( uint_t shell = 0; shell < numberOfShells( nRad, level, polynomialDegreeOfBasisFunctions< FunctionType >() ); shell++ )
   {
      VTKPointCloudOutput vtk( ".", "point_cloud_level_" + std::to_string( level ) + "_shell_" + std::to_string( shell ) );
      vtk.setPoints( shellData.points( shell ) );
      vtk.setValues( "u", shellData.values( "u", 0, shell ) );
      vtk.write();
   }

   const bool vtk = false;
   if ( vtk )
   {
      VTKOutput vtkMesh( ".", "data_with_mesh", storage );
      vtkMesh.add( u );
      vtkMesh.write( level );
   }
}

template < typename FunctionType >
void testRadialIntegerIDOutput( uint_t nTan, uint_t nRad, real_t rMin, real_t rMax, uint_t level )
{
   auto storage = setupSphericalShellStorage( nTan, nRad, rMin, rMax );

   FunctionType u( "u", storage, level, level );

   interpolateRadialShellID( u, rMin, rMax, nRad, level );

   VTKOutput vtkMesh( ".", "data_integer_id", storage );
   vtkMesh.add( u );
   vtkMesh.write( level );
}

template < typename FunctionType >
void testRadialIntegerIDOutput( uint_t nTan, std::vector<real_t> layers, uint_t level )
{
   auto storage = setupSphericalShellStorage( nTan, layers );

   FunctionType u( "u", storage, level, level );

   interpolateRadialShellID( u, layers, level );

   VTKOutput vtkMesh( ".", "variable_data_integer_id", storage );
   vtkMesh.add( u );
   vtkMesh.write( level );
}

} // namespace terraneo

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "testShellMath()" )
   terraneo::testShellMath();

   WALBERLA_LOG_INFO_ON_ROOT( "testRadialOutput< P1Function< real_t > >()" )
   terraneo::testRadialPointCloudOutput< P1Function< real_t > >( 5, 3, 0.5, 1.0, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "testRadialOutput< P1VectorFunction< real_t > >()" )
   terraneo::testRadialPointCloudOutput< P1VectorFunction< real_t > >( 5, 3, 0.5, 1.0, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "testRadialOutput< P2Function< real_t > >()" )
   terraneo::testRadialPointCloudOutput< P2Function< real_t > >( 5, 3, 0.5, 1.0, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "testRadialOutput< P2VectorFunction< real_t > >()" )
   terraneo::testRadialPointCloudOutput< P2VectorFunction< real_t > >( 5, 3, 0.5, 1.0, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "testRadialOutput< P1Function< real_t > >()" )
   terraneo::testRadialIntegerIDOutput< P1Function< int32_t > >( 5, 3, 0.5, 1.0, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "testRadialOutput< P2Function< real_t > >()" )
   terraneo::testRadialIntegerIDOutput< P2Function< int32_t > >( 5, 3, 0.5, 1.0, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "testVariableRadialOutput< P1Function< real_t > >()" )
   terraneo::testRadialIntegerIDOutput< P1Function< int32_t > >( 5, {0.5, 0.55, 1.0}, 2 );

   WALBERLA_LOG_INFO_ON_ROOT( "testVariableRadialOutput< P2Function< real_t > >()" )
   terraneo::testRadialIntegerIDOutput< P2Function< int32_t > >( 5, {0.5, 0.55, 1.0}, 2 );

   return 0;
}