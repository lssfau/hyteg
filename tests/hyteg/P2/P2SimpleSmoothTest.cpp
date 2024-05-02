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
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2MacroFace.hpp"
#include "hyteg/p2functionspace/P2MacroVertex.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constant_stencil_operator/P2ConstantOperator.hpp"
#include "constant_stencil_operator/P2MacroEdge.hpp"

namespace hyteg {

/// this test uses stencil weights which results in all values being one after the smoothing step so we can check easily
static void testP2Smooth()
{
   const uint_t level = 3;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../meshes/quad_2el.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto x   = std::make_shared< hyteg::P2Function< real_t > >( "x", storage, level, level );
   auto rhs = std::make_shared< hyteg::P2Function< real_t > >( "rhs", storage, level, level );

   P2ConstantLaplaceOperator p2operator( storage, level, level );

   typedef stencilDirection sD;

   real_t* vertexToVertexStencil;
   real_t* edgeToVertexStencil;

   real_t* edgeToEdgeStencil;
   real_t* vertexToEdgeStencil;

   for( auto& faceIt : storage->getFaces() )
   {
      Face& face = *faceIt.second;

      vertexToVertexStencil = face.getData( p2operator.getVertexToVertexOpr().getFaceStencilID() )->getPointer( level );
      edgeToVertexStencil   = face.getData( p2operator.getEdgeToVertexOpr().getFaceStencilID() )->getPointer( level );

      edgeToEdgeStencil   = face.getData( p2operator.getEdgeToEdgeOpr().getFaceStencilID() )->getPointer( level );
      vertexToEdgeStencil = face.getData( p2operator.getVertexToEdgeOpr().getFaceStencilID() )->getPointer( level );

      /// vertex dofs
      for( uint_t k = 0; k < vertexdof::macroface::neighborsWithCenter.size(); ++k )
      {
         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithCenter[k] )] = 1;
      }
      vertexToVertexStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] = -17.;

      for( uint_t k = 0; k < edgedof::macroface::neighborsFromVertex.size(); ++k )
      {
         edgeToVertexStencil[edgedof::stencilIndexFromVertex( edgedof::macroface::neighborsFromVertex[k] )] = 1;
      }

      /// horizontal edges
      for( uint_t k = 0; k < edgedof::macroface::neighborsFromHorizontalEdge.size(); ++k )
      {
         edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( edgedof::macroface::neighborsFromHorizontalEdge[k] )] = 1;
      }
      edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C )] = -7;

      for( uint_t k = 0; k < vertexdof::macroface::neighborsFromHorizontalEdge.size(); ++k )
      {
         vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( vertexdof::macroface::neighborsFromHorizontalEdge[k] )] =
             1;
      }

      /// diagonal edges
      for( uint_t k = 0; k < edgedof::macroface::neighborsFromDiagonalEdge.size(); ++k )
      {
         edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge( edgedof::macroface::neighborsFromDiagonalEdge[k] )] = 1;
      }
      edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge( stencilDirection::EDGE_DI_C )] = -7;

      for( uint_t k = 0; k < vertexdof::macroface::neighborsFromDiagonalEdge.size(); ++k )
      {
         vertexToEdgeStencil[vertexdof::stencilIndexFromDiagonalEdge( vertexdof::macroface::neighborsFromDiagonalEdge[k] )] = 1;
      }

      /// vertical edges
      for( uint_t k = 0; k < edgedof::macroface::neighborsFromVerticalEdge.size(); ++k )
      {
         edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge( edgedof::macroface::neighborsFromVerticalEdge[k] )] = 1;
      }
      edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge( stencilDirection::EDGE_VE_C )] = -7;

      for( uint_t k = 0; k < vertexdof::macroface::neighborsFromVerticalEdge.size(); ++k )
      {
         vertexToEdgeStencil[vertexdof::stencilIndexFromVerticalEdge( vertexdof::macroface::neighborsFromVerticalEdge[k] )] = 1;
      }

      std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return 1; };
      std::function< real_t( const hyteg::Point3D& ) > two  = []( const hyteg::Point3D& ) { return 2; };

      x->interpolate( ones, level );
      rhs->interpolate( ones, level );

      hyteg::communication::syncFunctionBetweenPrimitives( ( *x ), level );

      P2::macroface::smoothGaussSeidel( level,
                                        face,
                                        p2operator.getVertexToVertexOpr().getFaceStencilID(),
                                        p2operator.getEdgeToVertexOpr().getFaceStencilID(),
                                        x->getVertexDoFFunction().getFaceDataID(),
                                        p2operator.getVertexToEdgeOpr().getFaceStencilID(),
                                        p2operator.getEdgeToEdgeOpr().getFaceStencilID(),
                                        x->getEdgeDoFFunction().getFaceDataID(),
                                        rhs->getVertexDoFFunction().getFaceDataID(),
                                        rhs->getEdgeDoFFunction().getFaceDataID() );

      P2::macroface::smoothGaussSeidel( level,
                                        face,
                                        p2operator.getVertexToVertexOpr().getFaceStencilID(),
                                        p2operator.getEdgeToVertexOpr().getFaceStencilID(),
                                        x->getVertexDoFFunction().getFaceDataID(),
                                        p2operator.getVertexToEdgeOpr().getFaceStencilID(),
                                        p2operator.getEdgeToEdgeOpr().getFaceStencilID(),
                                        x->getEdgeDoFFunction().getFaceDataID(),
                                        rhs->getVertexDoFFunction().getFaceDataID(),
                                        rhs->getEdgeDoFFunction().getFaceDataID() );

      real_t* edgeDoFData   = face.getData( x->getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
      real_t* vertexDoFData = face.getData( x->getVertexDoFFunction().getFaceDataID() )->getPointer( level );

      for( const auto& it : hyteg::vertexdof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             vertexDoFData[hyteg::vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), sD::VERTEX_C )],
             1.,
             it.x() << " " << it.y() );
      }

      for( const auto& it : hyteg::edgedof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hyteg::edgedof::macroface::indexFromVertex( level, it.x(), it.y(), sD::EDGE_HO_E )],
             1.,
             it.x() << " " << it.y() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hyteg::edgedof::macroface::indexFromVertex( level, it.x(), it.y(), sD::EDGE_DI_NE )],
             1.,
             it.x() << " " << it.y() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hyteg::edgedof::macroface::indexFromVertex( level, it.x(), it.y(), sD::EDGE_VE_N )],
             1.,
             it.x() << " " << it.y() );
      }
   }

   /// check only the edge with two faces
   Edge* doubleEdge = nullptr;
   for( auto e : storage->getEdges() )
   {
      if( e.second->getNumNeighborFaces() == 2 )
      {
         doubleEdge = e.second.get();
      }
   }
   WALBERLA_CHECK( doubleEdge != nullptr );
   vertexToVertexStencil = doubleEdge->getData( p2operator.getVertexToVertexOpr().getEdgeStencilID() )->getPointer( level );
   edgeToVertexStencil   = doubleEdge->getData( p2operator.getEdgeToVertexOpr().getEdgeStencilID() )->getPointer( level );

   edgeToEdgeStencil   = doubleEdge->getData( p2operator.getEdgeToEdgeOpr().getEdgeStencilID() )->getPointer( level );
   vertexToEdgeStencil = doubleEdge->getData( p2operator.getVertexToEdgeOpr().getEdgeStencilID() )->getPointer( level );

   /// vertex dofs
   for( uint_t k = 0; k < 7; ++k )
   {
      vertexToVertexStencil[k] = 1;
   }
   vertexToVertexStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] = -17.;

   for( uint_t k = 0; k < 12; ++k )
   {
      edgeToVertexStencil[k] = 1;
   }

   /// edge dofs
   for( uint_t k = 0; k < 4; ++k )
   {
      vertexToEdgeStencil[k] = 1;
   }

   for( uint_t k = 0; k < 5; ++k )
   {
      edgeToEdgeStencil[k] = 1;
   }
   edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C )] = -7.;

   P2::macroedge::smoothGaussSeidel( level,
                                     *doubleEdge,
                                     p2operator.getVertexToVertexOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToVertexOpr().getEdgeStencilID(),
                                     x->getVertexDoFFunction().getEdgeDataID(),
                                     p2operator.getVertexToEdgeOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToEdgeOpr().getEdgeStencilID(),
                                     x->getEdgeDoFFunction().getEdgeDataID(),
                                     rhs->getVertexDoFFunction().getEdgeDataID(),
                                     rhs->getEdgeDoFFunction().getEdgeDataID() );

   P2::macroedge::smoothGaussSeidel( level,
                                     *doubleEdge,
                                     p2operator.getVertexToVertexOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToVertexOpr().getEdgeStencilID(),
                                     x->getVertexDoFFunction().getEdgeDataID(),
                                     p2operator.getVertexToEdgeOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToEdgeOpr().getEdgeStencilID(),
                                     x->getEdgeDoFFunction().getEdgeDataID(),
                                     rhs->getVertexDoFFunction().getEdgeDataID(),
                                     rhs->getEdgeDoFFunction().getEdgeDataID() );

   real_t* edgeDoFData   = doubleEdge->getData( x->getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );
   real_t* vertexDoFData = doubleEdge->getData( x->getVertexDoFFunction().getEdgeDataID() )->getPointer( level );

   for( const auto& it : hyteg::vertexdof::macroedge::Iterator( level, 0 ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( vertexDoFData[hyteg::vertexdof::macroedge::indexFromVertex( level, it.x(), sD::VERTEX_C )],
                                  1.,
                                  it.x() << " " << it.y() );
   }

   for( const auto& it : hyteg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFData[hyteg::edgedof::macroedge::indexFromVertex( level, it.x(), sD::EDGE_HO_E )],
                                  1.,
                                  it.x() << " " << it.y() );
   }
}

static void testP2JacobiSmooth()
{
   const uint_t level = 3;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../meshes/quad_2el.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto x   = std::make_shared< P2Function< real_t > >( "x", storage, level, level );
   auto tmp = std::make_shared< P2Function< real_t > >( "x", storage, level, level );
   auto rhs = std::make_shared< P2Function< real_t > >( "rhs", storage, level, level );

   P2ConstantLaplaceOperator p2operator( storage, level, level );

   typedef stencilDirection sD;

   real_t* vertexToVertexStencil;
   real_t* edgeToVertexStencil;

   real_t* edgeToEdgeStencil;
   real_t* vertexToEdgeStencil;

   std::function< real_t( const hyteg::Point3D& ) >                          ones = []( const hyteg::Point3D& ) { return 13; };
   std::function< real_t( const hyteg::Point3D& ) >                          two  = []( const hyteg::Point3D& ) { return 2; };
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > onesExtended =
       [&ones]( const hyteg::Point3D& xx, const std::vector< real_t >& ) { return ones( xx ); };

   tmp->interpolate( ones, level );
   rhs->interpolate( ones, level );

   hyteg::communication::syncFunctionBetweenPrimitives( ( *tmp ), level );
   hyteg::communication::syncFunctionBetweenPrimitives( ( *rhs ), level );

   for( auto e : storage->getEdges() )
   {
      Edge* edge = e.second.get();
      vertexdof::macroedge::interpolate( level, *edge, x->getVertexDoFFunction().getEdgeDataID(), {}, onesExtended );
      edgedof::macroedge::interpolate( level, *edge, x->getEdgeDoFFunction().getEdgeDataID(), {}, onesExtended );
   }
   for( auto v : storage->getVertices() )
   {
      Vertex* vertex = v.second.get();
      vertexdof::macrovertex::interpolate( *vertex, x->getVertexDoFFunction().getVertexDataID(), {}, onesExtended, level );
   }

   hyteg::communication::syncFunctionBetweenPrimitives( ( *x ), level );

   for( auto faceIt : storage->getFaces() )
   {
      Face* face = faceIt.second.get();

      vertexToVertexStencil = face->getData( p2operator.getVertexToVertexOpr().getFaceStencilID() )->getPointer( level );
      edgeToVertexStencil   = face->getData( p2operator.getEdgeToVertexOpr().getFaceStencilID() )->getPointer( level );

      edgeToEdgeStencil   = face->getData( p2operator.getEdgeToEdgeOpr().getFaceStencilID() )->getPointer( level );
      vertexToEdgeStencil = face->getData( p2operator.getVertexToEdgeOpr().getFaceStencilID() )->getPointer( level );

      /// vertex dofs
      for( uint_t k = 0; k < vertexdof::macroface::neighborsWithCenter.size(); ++k )
      {
         vertexToVertexStencil[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithCenter[k] )] = 1;
      }
      vertexToVertexStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] = -17.;

      for( uint_t k = 0; k < edgedof::macroface::neighborsFromVertex.size(); ++k )
      {
         edgeToVertexStencil[edgedof::stencilIndexFromVertex( edgedof::macroface::neighborsFromVertex[k] )] = 1;
      }

      /// horizontal edges
      for( uint_t k = 0; k < edgedof::macroface::neighborsFromHorizontalEdge.size(); ++k )
      {
         edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( edgedof::macroface::neighborsFromHorizontalEdge[k] )] = 1;
      }
      edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C )] = -7;

      for( uint_t k = 0; k < vertexdof::macroface::neighborsFromHorizontalEdge.size(); ++k )
      {
         vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge( vertexdof::macroface::neighborsFromHorizontalEdge[k] )] =
             1;
      }

      /// diagonal edges
      for( uint_t k = 0; k < edgedof::macroface::neighborsFromDiagonalEdge.size(); ++k )
      {
         edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge( edgedof::macroface::neighborsFromDiagonalEdge[k] )] = 1;
      }
      edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge( stencilDirection::EDGE_DI_C )] = -7;

      for( uint_t k = 0; k < vertexdof::macroface::neighborsFromDiagonalEdge.size(); ++k )
      {
         vertexToEdgeStencil[vertexdof::stencilIndexFromDiagonalEdge( vertexdof::macroface::neighborsFromDiagonalEdge[k] )] = 1;
      }

      /// vertical edges
      for( uint_t k = 0; k < edgedof::macroface::neighborsFromVerticalEdge.size(); ++k )
      {
         edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge( edgedof::macroface::neighborsFromVerticalEdge[k] )] = 1;
      }
      edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge( stencilDirection::EDGE_VE_C )] = -7;

      for( uint_t k = 0; k < vertexdof::macroface::neighborsFromVerticalEdge.size(); ++k )
      {
         vertexToEdgeStencil[vertexdof::stencilIndexFromVerticalEdge( vertexdof::macroface::neighborsFromVerticalEdge[k] )] = 1;
      }

      P2::macroface::smoothJacobiVertexDoF( level,
                                            *face,
                                            p2operator.getVertexToVertexOpr().getFaceStencilID(),
                                            tmp->getVertexDoFFunction().getFaceDataID(),
                                            x->getVertexDoFFunction().getFaceDataID(),
                                            p2operator.getEdgeToVertexOpr().getFaceStencilID(),
                                            tmp->getEdgeDoFFunction().getFaceDataID(),
                                            rhs->getVertexDoFFunction().getFaceDataID() );

      P2::macroface::smoothJacobiEdgeDoF( level,
                                          *face,
                                          p2operator.getVertexToEdgeOpr().getFaceStencilID(),
                                          tmp->getVertexDoFFunction().getFaceDataID(),
                                          p2operator.getEdgeToEdgeOpr().getFaceStencilID(),
                                          tmp->getEdgeDoFFunction().getFaceDataID(),
                                          x->getEdgeDoFFunction().getFaceDataID(),
                                          rhs->getEdgeDoFFunction().getFaceDataID() );

      real_t* edgeDoFData   = face->getData( x->getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
      real_t* vertexDoFData = face->getData( x->getVertexDoFFunction().getFaceDataID() )->getPointer( level );

      for( const auto& it : hyteg::vertexdof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             vertexDoFData[hyteg::vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), sD::VERTEX_C )],
             13.,
             it.x() << " " << it.y() );
      }

      for( const auto& it : hyteg::edgedof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hyteg::edgedof::macroface::indexFromVertex( level, it.x(), it.y(), sD::EDGE_HO_E )],
             13.,
             it.x() << " " << it.y() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hyteg::edgedof::macroface::indexFromVertex( level, it.x(), it.y(), sD::EDGE_DI_NE )],
             13.,
             it.x() << " " << it.y() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hyteg::edgedof::macroface::indexFromVertex( level, it.x(), it.y(), sD::EDGE_VE_N )],
             13.,
             it.x() << " " << it.y() );
      }
   }

   /// check only the edge with two faces
   Edge* doubleEdge = nullptr;
   for( auto e : storage->getEdges() )
   {
      if( e.second->getNumNeighborFaces() == 2 )
      {
         doubleEdge = e.second.get();
      }
   }
   WALBERLA_CHECK( doubleEdge != nullptr );
   vertexToVertexStencil = doubleEdge->getData( p2operator.getVertexToVertexOpr().getEdgeStencilID() )->getPointer( level );
   edgeToVertexStencil   = doubleEdge->getData( p2operator.getEdgeToVertexOpr().getEdgeStencilID() )->getPointer( level );

   edgeToEdgeStencil   = doubleEdge->getData( p2operator.getEdgeToEdgeOpr().getEdgeStencilID() )->getPointer( level );
   vertexToEdgeStencil = doubleEdge->getData( p2operator.getVertexToEdgeOpr().getEdgeStencilID() )->getPointer( level );

   /// vertex dofs
   for( uint_t k = 0; k < 7; ++k )
   {
      vertexToVertexStencil[k] = 1;
   }
   vertexToVertexStencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] = -17.;

   for( uint_t k = 0; k < 12; ++k )
   {
      edgeToVertexStencil[k] = 1;
   }

   /// edge dofs
   for( uint_t k = 0; k < 4; ++k )
   {
      vertexToEdgeStencil[k] = 1;
   }

   for( uint_t k = 0; k < 5; ++k )
   {
      edgeToEdgeStencil[k] = 1;
   }
   edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge( stencilDirection::EDGE_HO_C )] = -7.;


   P2::macroedge::smoothGaussSeidel( level,
                                     *doubleEdge,
                                     p2operator.getVertexToVertexOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToVertexOpr().getEdgeStencilID(),
                                     x->getVertexDoFFunction().getEdgeDataID(),
                                     p2operator.getVertexToEdgeOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToEdgeOpr().getEdgeStencilID(),
                                     x->getEdgeDoFFunction().getEdgeDataID(),
                                     rhs->getVertexDoFFunction().getEdgeDataID(),
                                     rhs->getEdgeDoFFunction().getEdgeDataID() );

   P2::macroedge::smoothGaussSeidel( level,
                                     *doubleEdge,
                                     p2operator.getVertexToVertexOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToVertexOpr().getEdgeStencilID(),
                                     x->getVertexDoFFunction().getEdgeDataID(),
                                     p2operator.getVertexToEdgeOpr().getEdgeStencilID(),
                                     p2operator.getEdgeToEdgeOpr().getEdgeStencilID(),
                                     x->getEdgeDoFFunction().getEdgeDataID(),
                                     rhs->getVertexDoFFunction().getEdgeDataID(),
                                     rhs->getEdgeDoFFunction().getEdgeDataID() );

   ///TODO: enable once jacobi on macroedges is implemented
#if 0
  real_t *edgeDoFData = doubleEdge->getData(x->getEdgeDoFFunction().getEdgeDataID()).getPointer(level);
  real_t *vertexDoFData = doubleEdge->getData(x->getVertexDoFFunction().getEdgeDataID()).getPointer(level);

  for (const auto &it : hyteg::vertexdof::macroedge::Iterator(level, 0)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFData[hyteg::vertexdof::macroedge::indexFromVertex(level,it.x(), sD::VERTEX_C)],
      1.,
      it.x() << " " << it.y());
  }

  for (const auto &it : hyteg::edgedof::macroedge::Iterator(level, 0)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      edgeDoFData[hyteg::edgedof::macroedge::indexFromVertex(level,it.x(), sD::EDGE_HO_E)],
      1.,
      it.x() << " " << it.y());
  }
#endif
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP2Smooth();
   hyteg::testP2JacobiSmooth();

   return EXIT_SUCCESS;
}
