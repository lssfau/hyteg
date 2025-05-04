/*
 * Copyright (c) 2024 Andreas Burkhart.
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
#include "hyteg/volumedofspace/MicroElementDoFTypification.hpp"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using namespace hyteg;

class VertexIdentifier
{
 public:
   VertexIdentifier( const Point3D pos, const MicroDoFType type, real_t tolerance = real_c( 1e-8 ) )
   : pos_( pos )
   , type_( type )
   , tolerance_( tolerance )
   {}

   bool testVertex( const Point3D testPos )
   {
      Point3D diff = pos_ - testPos;
      return diff.norm() < tolerance_;
   }

   bool testVertex( const VertexIdentifier& VertexIdent ) { return testVertex( VertexIdent.getPos() ); }

   MicroDoFType getType() const { return type_; }
   Point3D      getPos() const { return pos_; }

 private:
   Point3D      pos_;
   MicroDoFType type_;
   real_t       tolerance_;
};

class EdgeIdentifier
{
 public:
   EdgeIdentifier( const Point3D pos1, const Point3D pos2, const MicroDoFType type, real_t tolerance = real_c( 1e-8 ) )
   : pos1_( pos1 )
   , pos2_( pos2 )
   , diff1_( pos2 - pos1 )
   , type_( type )
   , tolerance_( tolerance )
   , dotProductTolerance_( std::pow( diff1_.norm(), 2 ) + tolerance )
   {}

   bool testVertex( const Point3D testPos )
   {
      Point3D diff2 = testPos - pos1_;
      // check if on line && between edge points
      real_t dot = diff1_.dot( diff2 );
      if ( ( diff2.cross( diff1_ ).norm() < tolerance_ ) && ( dot > -tolerance_ ) &&
           ( diff1_.dot( diff2 ) < dotProductTolerance_ ) )
      {
         return true;
      }
      return false;
   }

   bool testVertex( const VertexIdentifier& VertexIdent ) { return testVertex( VertexIdent.getPos() ); }

   bool testEdge( const Point3D testPos1, const Point3D testPos2 ) { return testVertex( testPos1 ) && testVertex( testPos2 ); }

   bool testEdge( const EdgeIdentifier& EdgeIdent ) { return testEdge( EdgeIdent.getPos1(), EdgeIdent.getPos2() ); }

   MicroDoFType getType() const { return type_; }
   Point3D      getPos1() const { return pos1_; }
   Point3D      getPos2() const { return pos2_; }

 private:
   Point3D      pos1_;
   Point3D      pos2_;
   Point3D      diff1_;
   MicroDoFType type_;
   real_t       tolerance_;
   real_t       dotProductTolerance_;
};

class FaceIdentifier
{
 public:
   FaceIdentifier( const Point3D      pos1,
                   const Point3D      pos2,
                   const Point3D      pos3,
                   const MicroDoFType type,
                   real_t             tolerance = real_c( 1e-8 ) )
   : pos1_( pos1 )
   , pos2_( pos2 )
   , pos3_( pos3 )
   , AB_( pos2_ - pos1_ )
   , BC_( pos3_ - pos2_ )
   , CA_( pos1_ - pos3_ )
   , normal_( AB_.cross( pos3_ - pos1_ ) )
   , type_( type )
   , tolerance_( tolerance )
   {
      normal_.normalize();
      dist_                     = normal_.dot( pos1_ );
      doubleSignedAreaTriangle_ = AB_.cross( pos3_ - pos1_ ).dot( normal_ );
   }

   bool testVertex( const Point3D testPos )
   {
      // check if point is on the correct plane
      bool onPlane = ( std::abs( normal_.dot( testPos ) - dist_ ) < tolerance_ );

      if ( onPlane )
      {
         // signed area of the subtriangles ABP, BCP, CAP
         real_t doubleSignedAreaSubTriangle1 = AB_.cross( testPos - pos1_ ).dot( normal_ );
         real_t doubleSignedAreaSubTriangle2 = BC_.cross( testPos - pos2_ ).dot( normal_ );
         real_t doubleSignedAreaSubTriangle3 = CA_.cross( testPos - pos3_ ).dot( normal_ );

         // barycentric coordinates
         real_t c1 = doubleSignedAreaSubTriangle2 / doubleSignedAreaTriangle_;
         real_t c2 = doubleSignedAreaSubTriangle3 / doubleSignedAreaTriangle_;
         real_t c3 = doubleSignedAreaSubTriangle1 / doubleSignedAreaTriangle_;

         // check if all barycentric coordinates are between 0 and 1 and that their sum is 1
         if ( ( c1 >= -tolerance_ ) && ( c1 <= real_c( 1.0 ) + tolerance_ ) && ( c2 >= -tolerance_ ) &&
              ( c2 <= real_c( 1.0 ) + tolerance_ ) && ( c3 >= -tolerance_ ) && ( c3 <= real_c( 1.0 ) + tolerance_ ) &&
              ( abs( c1 + c2 + c3 - real_c( 1.0 ) ) < tolerance_ ) )
            return true;
      }

      return false;
   }

   bool testVertex( const VertexIdentifier& VertexIdent ) { return testVertex( VertexIdent.getPos() ); }

   bool testEdge( const Point3D testPos1, const Point3D testPos2 ) { return testVertex( testPos1 ) && testVertex( testPos2 ); }

   bool testEdge( const EdgeIdentifier& EdgeIdent ) { return testEdge( EdgeIdent.getPos1(), EdgeIdent.getPos2() ); }

   bool testFace( const Point3D testPos1, const Point3D testPos2, const Point3D testPos3 )
   {
      return testVertex( testPos1 ) && testVertex( testPos2 ) && testVertex( testPos3 );
   }

   bool testFace( const FaceIdentifier& FaceIdent )
   {
      return testFace( FaceIdent.getPos1(), FaceIdent.getPos2(), FaceIdent.getPos3() );
   }

   MicroDoFType getType() const { return type_; }
   Point3D      getPos1() const { return pos1_; }
   Point3D      getPos2() const { return pos2_; }
   Point3D      getPos3() const { return pos3_; }

 private:
   Point3D pos1_;
   Point3D pos2_;
   Point3D pos3_;

   Point3D AB_;
   Point3D BC_;
   Point3D CA_;

   Point3D normal_;

   real_t dist_;
   real_t doubleSignedAreaTriangle_;

   MicroDoFType type_;
   real_t       tolerance_;
};

class FaceMicroDoFIdentifier
{
 public:
   FaceMicroDoFIdentifier( const std::shared_ptr< PrimitiveStorage >& storage,
                           const Face&                                face,
                           real_t                                     tolerance = real_c( 1e-8 ) )
   {
      // get macro vertices
      std::vector< hyteg::PrimitiveID > neighborVertices;

      face.getNeighborVertices( neighborVertices );

      for ( auto& vertexId : neighborVertices )
      {
         Vertex& vertex = *( storage->getVertex( vertexId ) );

         auto vertexCoords = vertex.getCoordinates();

         auto localVertexID = face.vertex_index( vertexId );

         MacroVertexIdentifiers_.emplace_back(
             Point3D( vertexCoords[0], vertexCoords[1], vertexCoords[2] ),
             static_cast< MicroDoFType >( static_cast< size_t >( MacroVertex0 ) + static_cast< size_t >( localVertexID ) ),
             tolerance );
      }

      // get macro edges
      std::vector< hyteg::PrimitiveID > neighborEdges;

      face.getNeighborEdges( neighborEdges );

      for ( auto& edgeId : neighborEdges )
      {
         Edge& edge = *( storage->getEdge( edgeId ) );

         auto edgeCoords = edge.getCoordinates();

         auto localEdgeID = face.edge_index( edgeId );

         MacroEdgeIdentifiers_.emplace_back(
             Point3D( edgeCoords[0][0], edgeCoords[0][1], edgeCoords[0][2] ),
             Point3D( edgeCoords[1][0], edgeCoords[1][1], edgeCoords[1][2] ),
             static_cast< MicroDoFType >( static_cast< size_t >( MacroEdge0 ) + static_cast< size_t >( localEdgeID ) ),
             tolerance );
      }
   }

   MicroDoFType getVertexMicroDoFType( const Point3D testPos )
   {
      MicroDoFType type = MacroInner;

      for ( auto& e : MacroEdgeIdentifiers_ )
      {
         if ( e.testVertex( testPos ) )
         {
            type = e.getType();
         }
      }

      for ( auto& v : MacroVertexIdentifiers_ )
      {
         if ( v.testVertex( testPos ) )
         {
            type = v.getType();
         }
      }

      return type;
   }

   MicroDoFType getVertexMicroDoFType( const VertexIdentifier& VertexIdent )
   {
      return getVertexMicroDoFType( VertexIdent.getPos() );
   }

   MicroDoFType getEdgeMicroDoFType( const Point3D testPos1, const Point3D testPos2 )
   {
      MicroDoFType type = MacroInner;

      for ( auto& e : MacroEdgeIdentifiers_ )
      {
         if ( e.testEdge( testPos1, testPos2 ) )
         {
            type = e.getType();
         }
      }

      return type;
   }

   MicroDoFType getEdgeMicroDoFType( const EdgeIdentifier& EdgeIdent )
   {
      return getEdgeMicroDoFType( EdgeIdent.getPos1(), EdgeIdent.getPos2() );
   }

   facedof::MicroFaceBorderType
       getMicroElementBorderType( const Point3D testPos0, const Point3D testPos1, const Point3D testPos2 )
   {
      facedof::MicroFaceBorderType type = facedof::MicroFaceBorderType::BordersNoEdge;

      for ( auto& e : MacroEdgeIdentifiers_ )
      {
         if ( e.testVertex( testPos0 ) )
         {
            type = type | facedof::getMicroFaceBorderTypeFromMicroDoFType( e.getType() );
         }
         if ( e.testVertex( testPos1 ) )
         {
            type = type | facedof::getMicroFaceBorderTypeFromMicroDoFType( e.getType() );
         }
         if ( e.testVertex( testPos2 ) )
         {
            type = type | facedof::getMicroFaceBorderTypeFromMicroDoFType( e.getType() );
         }
      }

      return type;
   }

   facedof::MicroFaceBorderType getMicroElementBorderType( FaceIdentifier& FaceIdent )
   {
      return getMicroElementBorderType( FaceIdent.getPos1(), FaceIdent.getPos2(), FaceIdent.getPos3() );
   }

 private:
   std::vector< VertexIdentifier > MacroVertexIdentifiers_;
   std::vector< EdgeIdentifier >   MacroEdgeIdentifiers_;
};

class CellMicroDoFIdentifier
{
 public:
   CellMicroDoFIdentifier( const std::shared_ptr< PrimitiveStorage >& storage,
                           const Cell&                                cell,
                           real_t                                     tolerance = real_c( 1e-8 ) )
   {
      // get macro vertices
      std::vector< hyteg::PrimitiveID > neighborVertices;

      cell.getNeighborVertices( neighborVertices );

      for ( auto& vertexId : neighborVertices )
      {
         Vertex& vertex = *( storage->getVertex( vertexId ) );

         auto vertexCoords = vertex.getCoordinates();

         auto localVertexID = cell.getLocalVertexID( vertexId );

         MacroVertexIdentifiers_.emplace_back(
             Point3D( vertexCoords[0], vertexCoords[1], vertexCoords[2] ),
             static_cast< MicroDoFType >( static_cast< size_t >( MacroVertex0 ) + static_cast< size_t >( localVertexID ) ),
             tolerance );

         // WALBERLA_LOG_INFO( "Vertex" << localVertexID << ": (" << vertexCoords[0] << ", " << vertexCoords[1] << ", "
         //                                     << vertexCoords[2] << ")" )
      }

      // get macro edges
      std::vector< hyteg::PrimitiveID > neighborEdges;

      cell.getNeighborEdges( neighborEdges );

      for ( auto& edgeId : neighborEdges )
      {
         Edge& edge = *( storage->getEdge( edgeId ) );

         auto edgeCoords = edge.getCoordinates();

         auto localEdgeID = cell.getLocalEdgeID( edgeId );

         MacroEdgeIdentifiers_.emplace_back(
             Point3D( edgeCoords[0][0], edgeCoords[0][1], edgeCoords[0][2] ),
             Point3D( edgeCoords[1][0], edgeCoords[1][1], edgeCoords[1][2] ),
             static_cast< MicroDoFType >( static_cast< size_t >( MacroEdge0 ) + static_cast< size_t >( localEdgeID ) ),
             tolerance );

         // WALBERLA_LOG_INFO( "Edge" << localEdgeID << ": (" << edgeCoords[0][0] << ", " << edgeCoords[0][1] << ", "
         //                                   << edgeCoords[0][2] << ") ; (" << edgeCoords[1][0] << ", " << edgeCoords[1][1] << ", "
         //                                   << edgeCoords[1][2] << ")" )
      }

      // get macro faces
      std::vector< hyteg::PrimitiveID > neighborFaces;

      cell.getNeighborFaces( neighborFaces );

      for ( auto& faceId : neighborFaces )
      {
         Face& face = *( storage->getFace( faceId ) );

         auto faceCoords = face.getCoordinates();

         auto localFaceID = cell.getLocalFaceID( faceId );

         MacroFaceIdentifiers_.emplace_back(
             Point3D( faceCoords[0][0], faceCoords[0][1], faceCoords[0][2] ),
             Point3D( faceCoords[1][0], faceCoords[1][1], faceCoords[1][2] ),
             Point3D( faceCoords[2][0], faceCoords[2][1], faceCoords[2][2] ),
             static_cast< MicroDoFType >( static_cast< size_t >( MacroFace0 ) + static_cast< size_t >( localFaceID ) ),
             tolerance );

         // WALBERLA_LOG_INFO( "Face" << localFaceID << ": (" << faceCoords[0][0] << ", " << faceCoords[0][1] << ", "
         //                                   << faceCoords[0][2] << ") ; (" << faceCoords[1][0] << ", " << faceCoords[1][1] << ", "
         //                                   << faceCoords[1][2] << ") ; (" << faceCoords[2][0] << ", " << faceCoords[2][1] << ", "
         //                                   << faceCoords[2][2] << ")" )
      }
   }

   MicroDoFType getVertexMicroDoFType( const Point3D testPos )
   {
      MicroDoFType type = MacroInner;

      for ( auto& f : MacroFaceIdentifiers_ )
      {
         if ( f.testVertex( testPos ) )
         {
            type = f.getType();
         }
      }

      for ( auto& e : MacroEdgeIdentifiers_ )
      {
         if ( e.testVertex( testPos ) )
         {
            type = e.getType();
         }
      }

      for ( auto& v : MacroVertexIdentifiers_ )
      {
         if ( v.testVertex( testPos ) )
         {
            type = v.getType();
         }
      }

      return type;
   }

   MicroDoFType getVertexMicroDoFType( const VertexIdentifier& VertexIdent )
   {
      return getVertexMicroDoFType( VertexIdent.getPos() );
   }

   MicroDoFType getEdgeMicroDoFType( const Point3D testPos1, const Point3D testPos2 )
   {
      MicroDoFType type = MacroInner;

      for ( auto& f : MacroFaceIdentifiers_ )
      {
         if ( f.testEdge( testPos1, testPos2 ) )
         {
            type = f.getType();
         }
      }

      for ( auto& e : MacroEdgeIdentifiers_ )
      {
         if ( e.testEdge( testPos1, testPos2 ) )
         {
            type = e.getType();
         }
      }

      return type;
   }

   MicroDoFType getEdgeMicroDoFType( const EdgeIdentifier& EdgeIdent )
   {
      return getEdgeMicroDoFType( EdgeIdent.getPos1(), EdgeIdent.getPos2() );
   }

   MicroDoFType getFaceMicroDoFType( const Point3D testPos1, const Point3D testPos2, const Point3D testPos3 )
   {
      MicroDoFType type = MacroInner;

      for ( auto& f : MacroFaceIdentifiers_ )
      {
         if ( f.testFace( testPos1, testPos2, testPos3 ) )
         {
            type = f.getType();
         }
      }

      return type;
   }

   MicroDoFType getFaceMicroDoFType( const FaceIdentifier& FaceIdent )
   {
      return getFaceMicroDoFType( FaceIdent.getPos1(), FaceIdent.getPos2(), FaceIdent.getPos3() );
   }

   celldof::MicroCellBorderType
       getMicroElementBorderType( const Point3D testPos0, const Point3D testPos1, const Point3D testPos2, const Point3D testPos3 )
   {
      celldof::MicroCellBorderType type = celldof::MicroCellBorderType::BordersNoFace;

      for ( auto& f : MacroFaceIdentifiers_ )
      {
         if ( f.testVertex( testPos0 ) )
         {
            type = type | celldof::getMicroCellDoFTypeFromMicroCellBorderType( f.getType() );
         }
         if ( f.testVertex( testPos1 ) )
         {
            type = type | celldof::getMicroCellDoFTypeFromMicroCellBorderType( f.getType() );
         }
         if ( f.testVertex( testPos2 ) )
         {
            type = type | celldof::getMicroCellDoFTypeFromMicroCellBorderType( f.getType() );
         }
         if ( f.testVertex( testPos3 ) )
         {
            type = type | celldof::getMicroCellDoFTypeFromMicroCellBorderType( f.getType() );
         }
      }

      return type;
   }

 private:
   std::vector< VertexIdentifier > MacroVertexIdentifiers_;
   std::vector< EdgeIdentifier >   MacroEdgeIdentifiers_;
   std::vector< FaceIdentifier >   MacroFaceIdentifiers_;
};

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   if ( sizeof( real_t ) < 8 )
   {
      // The tolerances in this test are designed for at least double precision.
      WALBERLA_LOG_INFO_ON_ROOT( "Single precision or lower detected. Aborting test." );

      return EXIT_SUCCESS;
   }

   const uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   //const uint_t nProc = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   const bool Test2D = true;
   const bool Test3D = true;

   if ( Test2D )
   {
      uint_t minLevel_ = 0;
      for ( uint_t maxLevel_ = minLevel_; maxLevel_ <= 4; maxLevel_++ )
      {
         uint_t n_ = 1;

         // MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );
         MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISSCROSS, n_, n_ );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // set load balancing
         loadbalancing::roundRobinVolume( setupStorage );

         // set boundary flags
         setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

         auto storage_ = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

         BoundaryCondition bc = BoundaryCondition::create0123BC();

         // loop over all faces
         for ( auto& it : storage_->getFaces() )
         {
            Face& face = *it.second;

            // get macro face coordinates
            auto coords = face.getCoordinates();

            WALBERLA_LOG_INFO( "Rank " << rank << " Macro Face: "
                                       << "(" << coords[0][0] << ", " << coords[0][1] << ", " << coords[0][2] << ")"
                                       << " ; "
                                       << "(" << coords[1][0] << ", " << coords[1][1] << ", " << coords[1][2] << ")"
                                       << " ; "
                                       << "(" << coords[2][0] << ", " << coords[2][1] << ", " << coords[2][2] << ")" );

            FaceMicroDoFIdentifier MicroDoFIdentifier( storage_, face, real_c( 1e-8 ) );

            // create MicroDoFType to DoFType map
            auto MicroDoFTypeToDoFTypeMap = facedof::getMicroDoFTypeToDoFTypeMapFromMacroFace( storage_, face, bc );

            WALBERLA_LOG_INFO( "MicroDoFTypeToDoFTypeMap:" )
            for ( const auto& e : MicroDoFTypeToDoFTypeMap )
            {
               WALBERLA_LOG_INFO( "( " << e.first << ", " << e.second << " )" );
            }

            // loop over micro faces
            for ( const auto& faceType : facedof::allFaceTypes )
            {
               switch ( faceType )
               {
               case facedof::FaceType::BLUE:
                  WALBERLA_LOG_INFO( "FaceType BLUE:" );
                  break;

               case facedof::FaceType::GRAY:
                  WALBERLA_LOG_INFO( "FaceType GRAY:" );
                  break;
               }

               for ( const auto& microFace : facedof::macroface::Iterator( maxLevel_, faceType, 0 ) )
               {
                  std::array< indexing::Index, 3 > verts =
                      facedof::macroface::getMicroVerticesFromMicroFace( microFace, faceType );

                  std::array< Point3D, 3 > coords_micro;
                  for ( uint_t k = 0; k < 3; k++ )
                  {
                     coords_micro[k] = vertexdof::macroface::coordinateFromIndex( maxLevel_, face, verts[k] );
                  }

                  Point3D pos0( coords_micro[0][0], coords_micro[0][1], coords_micro[0][2] );
                  Point3D pos1( coords_micro[1][0], coords_micro[1][1], coords_micro[1][2] );
                  Point3D pos2( coords_micro[2][0], coords_micro[2][1], coords_micro[2][2] );

                  std::vector< MicroDoFType > empiricalMicroDoFTypes = { MicroDoFIdentifier.getVertexMicroDoFType( pos0 ),
                                                                         MicroDoFIdentifier.getVertexMicroDoFType( pos1 ),
                                                                         MicroDoFIdentifier.getVertexMicroDoFType( pos2 ),
                                                                         MicroDoFIdentifier.getEdgeMicroDoFType( pos0, pos1 ),
                                                                         MicroDoFIdentifier.getEdgeMicroDoFType( pos0, pos2 ),
                                                                         MicroDoFIdentifier.getEdgeMicroDoFType( pos1, pos2 ) };

                  auto MicroDoFTypes = getMicroDoFTypesFromMicroFace( faceType, microFace, maxLevel_ );

                  facedof::MicroFaceBorderType empiricalBorderType =
                      MicroDoFIdentifier.getMicroElementBorderType( pos0, pos1, pos2 );
                  WALBERLA_LOG_INFO( "empirical BorderType = " << empiricalBorderType );

                  std::vector< DoFType > EmpiricalFenicsOrderDofTypes =
                      facedof::getP2DoFTypesFromMicroDoFTypesFEniCSOrdering2D( MicroDoFTypeToDoFTypeMap, empiricalMicroDoFTypes );

                  std::vector< DoFType > FenicsOrderDofTypes =
                      facedof::getP2DoFTypesFromMicroDoFTypesFEniCSOrdering2D( MicroDoFTypeToDoFTypeMap, MicroDoFTypes );

                  WALBERLA_LOG_INFO(
                      "Rank " << rank << " Micro, Index "
                              << "(" << microFace[0] << ", " << microFace[1] << ", " << microFace[2] << ")"
                              << ": "
                              << "(" << coords_micro[0][0] << ", " << coords_micro[0][1] << ", " << coords_micro[0][2] << ")"
                              << " ; "
                              << "(" << coords_micro[1][0] << ", " << coords_micro[1][1] << ", " << coords_micro[1][2] << ")"
                              << " ; "
                              << "(" << coords_micro[2][0] << ", " << coords_micro[2][1] << ", " << coords_micro[2][2] << ")" );

                  std::stringstream EmpiricalMicroDoFString;
                  EmpiricalMicroDoFString << empiricalMicroDoFTypes[0];
                  for ( uint_t i = 1; i < empiricalMicroDoFTypes.size(); i++ )
                  {
                     EmpiricalMicroDoFString << ", " << empiricalMicroDoFTypes[i];
                  }
                  WALBERLA_LOG_INFO( "empiricalMicroDoFTypes " << EmpiricalMicroDoFString.str() );

                  std::stringstream MicroDoFString;
                  MicroDoFString << MicroDoFTypes[0];
                  for ( uint_t i = 1; i < MicroDoFTypes.size(); i++ )
                  {
                     MicroDoFString << ", " << MicroDoFTypes[i];
                  }
                  WALBERLA_LOG_INFO( "MicroDoFTypes: " << MicroDoFString.str() );

                  for ( uint_t i = 0; i < empiricalMicroDoFTypes.size(); i++ )
                  {
                     if ( empiricalMicroDoFTypes[i] != MicroDoFTypes[i] )
                     {
                        WALBERLA_ABORT( "MicroDoFTypeMismatch between index based and empirically determined MicroDoFTypes: "
                                        << empiricalMicroDoFTypes[i] << " != " << MicroDoFTypes[i] );
                     }
                  }

                  std::stringstream EmpiricalFenicsOrderString;
                  EmpiricalFenicsOrderString << EmpiricalFenicsOrderDofTypes[0];
                  for ( uint_t i = 1; i < EmpiricalFenicsOrderDofTypes.size(); i++ )
                  {
                     EmpiricalFenicsOrderString << ", " << EmpiricalFenicsOrderDofTypes[i];
                  }
                  WALBERLA_LOG_INFO( "Empirical Fenics ordering: " << EmpiricalFenicsOrderString.str() );

                  std::stringstream FenicsOrderString;
                  FenicsOrderString << FenicsOrderDofTypes[0];
                  for ( uint_t i = 1; i < FenicsOrderDofTypes.size(); i++ )
                  {
                     FenicsOrderString << ", " << FenicsOrderDofTypes[i];
                  }
                  WALBERLA_LOG_INFO( "Fenics ordering: " << FenicsOrderString.str() );

                  for ( uint_t i = 0; i < EmpiricalFenicsOrderDofTypes.size(); i++ )
                  {
                     if ( EmpiricalFenicsOrderDofTypes[i] != FenicsOrderDofTypes[i] )
                     {
                        WALBERLA_ABORT( "DoFTypeMismatch between index based and empirically determined DoFTypes: "
                                        << FenicsOrderDofTypes[i] << " != " << EmpiricalFenicsOrderDofTypes[i] );
                     }
                  }
               }
            }
         }
      }
   }

   if ( Test3D )
   {
      uint_t minLevel_ = 0;
      for ( uint_t maxLevel_ = minLevel_; maxLevel_ <= 4; maxLevel_++ )
      {
         uint_t n_ = 1;

         // MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
         MeshInfo              meshInfo = MeshInfo::meshCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), n_, n_, n_ );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // set load balancing
         loadbalancing::roundRobinVolume( setupStorage );

         // set boundary flags
         setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

         auto storage_ = std::make_shared< PrimitiveStorage >( setupStorage, 0 );

         BoundaryCondition bc = BoundaryCondition::create0123BC();

         // loop over all cells
         for ( auto& it : storage_->getCells() )
         {
            Cell& cell = *it.second;

            // get macro face coordinates
            auto coords = cell.getCoordinates();

            WALBERLA_LOG_INFO( "Rank " << rank << " Macro Face: "
                                       << "(" << coords[0][0] << ", " << coords[0][1] << ", " << coords[0][2] << ")"
                                       << " ; "
                                       << "(" << coords[1][0] << ", " << coords[1][1] << ", " << coords[1][2] << ")"
                                       << " ; "
                                       << "(" << coords[2][0] << ", " << coords[2][1] << ", " << coords[2][2] << ")"
                                       << " ; "
                                       << "(" << coords[3][0] << ", " << coords[3][1] << ", " << coords[3][2] << ")" );

            CellMicroDoFIdentifier MicroDoFIdentifier( storage_, cell, real_c( 1e-8 ) );

            // create MicroDoFType to DoFType map
            auto MicroDoFTypeToDoFTypeMap = celldof::getMicroDoFTypeToDoFTypeMapFromMacroCell( storage_, cell, bc );

            // loop over micro cells
            for ( const auto& cellType : celldof::allCellTypes )
            {
               WALBERLA_LOG_INFO( "CellType: " << celldof::CellTypeToStr.at( cellType ) );
               for ( const auto& microCell : celldof::macrocell::Iterator( maxLevel_, cellType, 0 ) )
               {
                  std::array< indexing::Index, 4 > verts =
                      celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cellType );
                  std::array< Point3D, 4 > coords_micro;
                  for ( uint_t k = 0; k < 4; ++k )
                  {
                     coords_micro[k] = vertexdof::macrocell::coordinateFromIndex( maxLevel_, cell, verts[k] );
                  }

                  WALBERLA_LOG_INFO(
                      "Rank " << rank << " Micro, Index "
                              << "(" << microCell[0] << ", " << microCell[1] << ", " << microCell[2] << ")"
                              << ": "
                              << "(" << coords_micro[0][0] << ", " << coords_micro[0][1] << ", " << coords_micro[0][2] << ")"
                              << " ; "
                              << "(" << coords_micro[1][0] << ", " << coords_micro[1][1] << ", " << coords_micro[1][2] << ")"
                              << " ; "
                              << "(" << coords_micro[2][0] << ", " << coords_micro[2][1] << ", " << coords_micro[2][2] << ")"
                              << " ; "
                              << "(" << coords_micro[3][0] << ", " << coords_micro[3][1] << ", " << coords_micro[3][2] << ")" );

                  Point3D pos0( coords_micro[0][0], coords_micro[0][1], coords_micro[0][2] );
                  Point3D pos1( coords_micro[1][0], coords_micro[1][1], coords_micro[1][2] );
                  Point3D pos2( coords_micro[2][0], coords_micro[2][1], coords_micro[2][2] );
                  Point3D pos3( coords_micro[3][0], coords_micro[3][1], coords_micro[3][2] );

                  auto empiricalBorderType = MicroDoFIdentifier.getMicroElementBorderType( pos0, pos1, pos2, pos3 );
                  WALBERLA_LOG_INFO( "empirical BorderType = " << empiricalBorderType );

                  auto MicroDoFTypes = getMicroDoFTypesFromMicroCell( cellType, microCell, maxLevel_ );

                  std::vector< MicroDoFType > empiricalMicroDoFTypes = {
                      MicroDoFIdentifier.getVertexMicroDoFType( pos0 ),
                      MicroDoFIdentifier.getVertexMicroDoFType( pos1 ),
                      MicroDoFIdentifier.getVertexMicroDoFType( pos2 ),
                      MicroDoFIdentifier.getVertexMicroDoFType( pos3 ),
                      MicroDoFIdentifier.getEdgeMicroDoFType( pos0, pos1 ),
                      MicroDoFIdentifier.getEdgeMicroDoFType( pos0, pos2 ),
                      MicroDoFIdentifier.getEdgeMicroDoFType( pos1, pos2 ),
                      MicroDoFIdentifier.getEdgeMicroDoFType( pos0, pos3 ),
                      MicroDoFIdentifier.getEdgeMicroDoFType( pos1, pos3 ),
                      MicroDoFIdentifier.getEdgeMicroDoFType( pos2, pos3 ),
                      MicroDoFIdentifier.getFaceMicroDoFType( pos0, pos1, pos2 ),
                      MicroDoFIdentifier.getFaceMicroDoFType( pos0, pos1, pos3 ),
                      MicroDoFIdentifier.getFaceMicroDoFType( pos0, pos2, pos3 ),
                      MicroDoFIdentifier.getFaceMicroDoFType( pos1, pos2, pos3 ) };

                  std::stringstream EmpiricalMicroDoFString;
                  EmpiricalMicroDoFString << empiricalMicroDoFTypes[0];
                  for ( uint_t i = 1; i < empiricalMicroDoFTypes.size(); i++ )
                  {
                     EmpiricalMicroDoFString << ", " << empiricalMicroDoFTypes[i];
                  }
                  WALBERLA_LOG_INFO( "empiricalMicroDoFTypes " << EmpiricalMicroDoFString.str() );

                  std::stringstream MicroDoFString;
                  MicroDoFString << MicroDoFTypes[0];
                  for ( uint_t i = 1; i < MicroDoFTypes.size(); i++ )
                  {
                     MicroDoFString << ", " << MicroDoFTypes[i];
                  }
                  WALBERLA_LOG_INFO( "MicroDoFTypes " << MicroDoFString.str() );

                  for ( uint_t i = 0; i < empiricalMicroDoFTypes.size(); i++ )
                  {
                     if ( empiricalMicroDoFTypes[i] != MicroDoFTypes[i] )
                     {
                        WALBERLA_ABORT( "MicroDoFTypeMismatch between index based and empirically determined MicroDoFTypes: "
                                        << empiricalMicroDoFTypes[i] << " != " << MicroDoFTypes[i] );
                     }
                  }

                  std::vector< DoFType > EmpiricalFenicsOrderDofTypes =
                      celldof::getP2DoFTypesFromMicroDoFTypesFEniCSOrdering3D( MicroDoFTypeToDoFTypeMap, empiricalMicroDoFTypes );

                  std::vector< DoFType > FenicsOrderDofTypes =
                      celldof::getP2DoFTypesFromMicroDoFTypesFEniCSOrdering3D( MicroDoFTypeToDoFTypeMap, MicroDoFTypes );

                  std::stringstream EmpiricalFenicsOrderString;
                  EmpiricalFenicsOrderString << EmpiricalFenicsOrderDofTypes[0];
                  for ( uint_t i = 1; i < EmpiricalFenicsOrderDofTypes.size(); i++ )
                  {
                     EmpiricalFenicsOrderString << ", " << EmpiricalFenicsOrderDofTypes[i];
                  }
                  WALBERLA_LOG_INFO( "Empirical Fenics ordering: " << EmpiricalFenicsOrderString.str() );

                  std::stringstream FenicsOrderString;
                  FenicsOrderString << FenicsOrderDofTypes[0];
                  for ( uint_t i = 1; i < FenicsOrderDofTypes.size(); i++ )
                  {
                     FenicsOrderString << ", " << FenicsOrderDofTypes[i];
                  }
                  WALBERLA_LOG_INFO( "Fenics ordering: " << FenicsOrderString.str() );

                  for ( uint_t i = 0; i < EmpiricalFenicsOrderDofTypes.size(); i++ )
                  {
                     if ( EmpiricalFenicsOrderDofTypes[i] != FenicsOrderDofTypes[i] )
                     {
                        WALBERLA_ABORT( "DoFTypeMismatch between index based and empirically determined DoFTypes: "
                                        << FenicsOrderDofTypes[i] << " != " << EmpiricalFenicsOrderDofTypes[i] );
                     }
                  }
               }
            }
         }
      }
   }

   return EXIT_SUCCESS;
}