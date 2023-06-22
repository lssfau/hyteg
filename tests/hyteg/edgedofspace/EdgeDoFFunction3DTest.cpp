/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include "core/debug/all.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t level = 2;

   EdgeDoFFunction< real_t > x1( "x1", storage, level, level );
   EdgeDoFFunction< real_t > x2( "x2", storage, level, level );
   EdgeDoFFunction< real_t > xSum( "xSum", storage, level, level );

   std::function< real_t( const Point3D& ) > fn1 = []( const Point3D& xx ) { return 2 * xx[0] + xx[1] * xx[2]; };

   std::function< real_t( const Point3D& ) > fn2 = []( const Point3D& xx ) { return xx[0] + 2 * xx[1] * xx[2]; };

   std::function< real_t( const Point3D& ) > fnSum = []( const Point3D& xx ) { return 3 * xx[0] + 3 * xx[1] * xx[2]; };

   std::function< real_t( const Point3D& ) > zero = []( const Point3D& ) { return real_c( 0 ); };

   std::function< real_t( const Point3D& ) > ones = []( const Point3D& ) { return real_c( 1 ); };

//   VTKOutput vtkOutput( "../../output", "EdgeDoFInterpolation3DTest", storage );
//   vtkOutput.add( x1 );
//   vtkOutput.add( x2 );
//   vtkOutput.add( xSum );
//
//   vtkOutput.write( level, 0 );

   x1.interpolate( fn1, level );
   xSum.interpolate( fn2, level );
   x2.assign( {1.0}, {xSum}, level );
   xSum.interpolate( zero, level );
   xSum.add( {1.0, 1.0}, {x1, x2}, level );
   x1.interpolate( fnSum, level );

   x2.interpolate( ones, level );
   const real_t dotProduct = x2.dotLocal( x2, level );

//   vtkOutput.write( level, 1 );

   for ( const auto& itCell : storage->getCells() )
   {
      const auto sumData  = itCell.second->getData( xSum.getCellDataID() )->getPointer( level );
      const auto testData = itCell.second->getData( x1.getCellDataID() )->getPointer( level );

      const auto length = itCell.second->getData( xSum.getCellDataID() )->getSize( level );

      for ( uint_t i = 0; i < length; i++ )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( sumData[i], testData[i] );
      }
   }

   WALBERLA_CHECK_FLOAT_EQUAL( dotProduct, 6 * real_c( 20 ) + real_c( 10 ) );

   x1.interpolate( real_c( 1.0 / 3.0 ), level, All );
   x1.invertElementwise( level, All );
   WALBERLA_CHECK_FLOAT_EQUAL( x1.getMaxValue( level, All ), real_c( 3 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( x1.getMinValue( level, All ), real_c( 3 ) );
}
