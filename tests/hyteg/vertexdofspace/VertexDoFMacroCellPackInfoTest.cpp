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
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"

namespace hyteg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;

static void testVertexDoFMacroCellPackInfo( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode )
{
  const uint_t level = 4;

  auto storage = PrimitiveStorage::createFromGmshFile( "../../meshes/3D/tet_1el.msh" );

  auto x = std::make_shared< vertexdof::VertexDoFFunction< real_t > >( "x", storage, level, level );

  // Writing 1.0 to all macro-faces
  for ( const auto & f : storage->getFaces() )
  {
    auto faceData = f.second->getData( x->getFaceDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macroface::Iterator( level ) )
    {
      faceData[vertexdof::macroface::indexFromVertex( level, it.x(), it.y(),
                                                      stencilDirection::VERTEX_C )] = 1.0;
    }
  }

  // Macro-cell should still be zero everywhere - particularly also at the borders
  for ( const auto & f : storage->getCells() )
  {
    auto cellData = f.second->getData( x->getCellDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macrocell::BoundaryIterator( level, 0, 1, 2 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )], 0.0 );
    }
    for ( const auto & it : vertexdof::macrocell::BoundaryIterator( level, 1, 2, 3 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )], 0.0 );
    }
  }

  // Communicating macro-face data to the macro-cell
  x->setLocalCommunicationMode( localCommunicationMode );
  x->communicate< Face, Cell >( level );

  // Macro-cell DoFs are 1.0 at the borders now
  for ( const auto & f : storage->getCells() )
  {
    auto cellData = f.second->getData( x->getCellDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macrocell::BoundaryIterator( level, 0, 1, 2 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )], 1.0 );
    }
    for ( const auto & it : vertexdof::macrocell::BoundaryIterator( level, 1, 2, 3 ) )
    {
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C )], 1.0 );
    }
  }

  // Now we check the correct rotation during the unpacking process by interpolation of a function
  std::function< real_t( const hyteg::Point3D & ) > expr = []( const hyteg::Point3D & xx ) -> real_t { return real_c( (1.0L/2.0L)*sin(2*xx[0])*sinh(xx[1]) ) * real_c( xx[2] ); };
  x->interpolate( expr, level );

  x->communicate< Vertex, Edge >( level );
  x->communicate< Edge, Face >( level );
  x->communicate< Face, Cell >( level );

  for ( const auto & f : storage->getCells() )
  {
    auto cellData = f.second->getData( x->getCellDataID() )->getPointer( level );
    for ( const auto & it : vertexdof::macrocell::Iterator( level ) )
    {
      const Point3D coordinate = vertexdof::macrocell::coordinateFromIndex( level, *f.second, it );
      const uint_t  idx        = vertexdof::macrocell::indexFromVertex( level, it.x(), it.y(), it.z(), stencilDirection::VERTEX_C );
      WALBERLA_CHECK_FLOAT_EQUAL( cellData[ idx ], expr( coordinate ) );
    }
  }
}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testVertexDoFMacroCellPackInfo( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::BUFFERED_MPI );
   hyteg::testVertexDoFMacroCellPackInfo( hyteg::communication::BufferedCommunicator::LocalCommunicationMode::DIRECT );

   return EXIT_SUCCESS;
}
