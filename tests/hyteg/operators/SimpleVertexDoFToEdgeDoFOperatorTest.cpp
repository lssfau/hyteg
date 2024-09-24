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

#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_t;

static void testOperator()
{
   const uint_t level = 3;

   MeshInfo mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/quad_4el.msh" ) );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto vertexDof = std::make_shared< P1Function< real_t > >( "vertexDof", storage, level, level );
   auto edgeDof   = std::make_shared< EdgeDoFFunction< real_t > >( "edgeDof", storage, level, level );

   GenericVertexDoFToEdgeDoFOperator vertexToEdgeOperator( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1; };
   std::function< real_t( const hyteg::Point3D&, const std::vector< real_t >& ) > onesVec =
       []( const hyteg::Point3D&, const std::vector< real_t >& ) { return 1; };

   for ( auto& faceIT : storage->getFaces() )
   {
      auto         face            = faceIT.second;
      const uint_t faceStencilSize = face->getData( vertexToEdgeOperator.getFaceStencilID() )->getSize( level );
      real_t*      faceStencilData = face->getData( vertexToEdgeOperator.getFaceStencilID() )->getPointer( level );
      for ( uint_t i = 0; i < faceStencilSize; ++i )
      {
         faceStencilData[i] = 1;
      }
      hyteg::vertexdof::macroface::interpolate< real_t >( storage, level, *face, vertexDof->getFaceDataID(), {}, onesVec );
   }

   for ( auto& edgeIT : storage->getEdges() )
   {
      auto         edge            = edgeIT.second;
      const uint_t edgeStencilSize = edge->getData( vertexToEdgeOperator.getEdgeStencilID() )->getSize( level );
      real_t*      edgeStencilData = edge->getData( vertexToEdgeOperator.getEdgeStencilID() )->getPointer( level );
      for ( uint_t i = 0; i < edgeStencilSize; ++i )
      {
         edgeStencilData[i] = 1;
      }
      hyteg::vertexdof::macroedge::interpolate< real_t >( storage, level, *edge, vertexDof->getEdgeDataID(), {}, onesVec );
   }

   for ( auto& vertexIT : storage->getVertices() )
   {
      auto vertex = vertexIT.second;
      hyteg::vertexdof::macrovertex::interpolate< real_t >( storage, *vertex, vertexDof->getVertexDataID(), {}, onesVec, level );
   }

   vertexToEdgeOperator.apply( *vertexDof, *edgeDof, level, hyteg::All );

   // Pull all halos
   edgeDof->communicate< Edge, Face >( level );
   edgeDof->communicate< Face, Edge >( level );
   edgeDof->communicate< Edge, Vertex >( level );
   edgeDof->communicate< Vertex, Edge >( level );

   //  for (auto &faceIT : storage->getFaces()) {
   //    auto face = faceIT.second;
   //    hyteg::edgedof::macroface::printFunctionMemory<real_t, level>(*face, edgeDof->getFaceDataID());
   //  }
   //
   //  for (auto &edgeIT : storage->getEdges()) {
   //    auto edge = edgeIT.second;
   //    hyteg::edgedof::macroedge::printFunctionMemory<real_t, level>(*edge, edgeDof->getEdgeDataID());
   //  }

   for ( auto& faceIT : storage->getFaces() )
   {
      auto         face = faceIT.second;
      const uint_t size = face->getData( edgeDof->getFaceDataID() )->getSize( level );
      real_t*      data = face->getData( edgeDof->getFaceDataID() )->getPointer( level );
      for ( uint_t i = 0; i < size; ++i )
      {
         ///the values on boundary can be 4 or 3 depending wether there is an adjacent face or not
         ///this check could be better but would be much more complicated
         WALBERLA_CHECK( walberla::floatIsEqual( data[i], real_c( 4.0 ) ) || walberla::floatIsEqual( data[i], real_c( 3.0 ) ) );
      }
   }

   for ( auto& edgeIT : storage->getEdges() )
   {
      auto         edge = edgeIT.second;
      const uint_t size = edge->getData( edgeDof->getEdgeDataID() )->getSize( level );
      real_t*      data = edge->getData( edgeDof->getEdgeDataID() )->getPointer( level );
      for ( uint_t i = 0; i < size; ++i )
      {
         ///this check could also be imporoved
         WALBERLA_CHECK( walberla::floatIsEqual( data[i], real_c( 4.0 ) ) || walberla::floatIsEqual( data[i], real_c( 3.0 ) ) );
      }
      hyteg::vertexdof::macroedge::interpolate< real_t >( storage, level, *edge, vertexDof->getEdgeDataID(), {}, onesVec );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testOperator();

   return EXIT_SUCCESS;
}
