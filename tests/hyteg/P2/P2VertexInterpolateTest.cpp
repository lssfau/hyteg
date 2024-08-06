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
#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = prependHyTeGMeshDir( "tri_1el_neumann.msh" );

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   size_t level = 2;

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P2Function< real_t >                      u( "u", storage, level, level );
   std::function< real_t( const hyteg::Point3D& ) > testExpression = []( const hyteg::Point3D& x ) { return x[1] + 1.0; };
   u.interpolate( testExpression, level, hyteg::All );

   // Sync interpolated function values
   u.getEdgeDoFFunction().communicate< Vertex, Edge >( level );
   u.getEdgeDoFFunction().communicate< Edge, Face >( level );
   u.getEdgeDoFFunction().communicate< Face, Edge >( level );
   u.getEdgeDoFFunction().communicate< Edge, Vertex >( level );

   // We assume that the bottom left vertex has following connectivity
   //
   // EdgeDoF(1)
   //    |
   //    |   EdgeDoF(2)
   //    |   /
   // Vertex(0) --- EdgeDoF(0)

   // Get bottom left vertex
   auto vertex = storage->getVertex( PrimitiveID::create( 0 ) );

   // Get vertex values
   auto vertexEdgeData = vertex->getData( u.getEdgeDoFFunction().getVertexDataID() )->getPointer( level );

   WALBERLA_CHECK_FLOAT_EQUAL( vertexEdgeData[0], 1.0 );
   WALBERLA_CHECK_FLOAT_EQUAL( vertexEdgeData[1], 1.125 );
   WALBERLA_CHECK_FLOAT_EQUAL( vertexEdgeData[2], 1.125 );

   return EXIT_SUCCESS;
}
