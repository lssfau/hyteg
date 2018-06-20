#include <tinyhhg_core/communication/Syncing.hpp>

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"

#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hhg {

/// this test uses stencil weights which results in all values being one after the smoothing step so we can check easily
static void testP2Smooth()
{
   const uint_t level = 3;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/quad_2el.msh" );

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto x   = std::make_shared< P2Function< real_t > >( "x", storage, level, level );
   auto rhs = std::make_shared< P2Function< real_t > >( "rhs", storage, level, level );

   P2ConstantLaplaceOperator p2operator( storage, level, level );

   typedef stencilDirection sD;

   real_t* vertexToVertexStencil;
   real_t* edgeToVertexStencil;

   real_t* edgeToEdgeStencil;
   real_t* vertexToEdgeStencil;

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

      std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1; };
      std::function< real_t( const hhg::Point3D& ) > two  = []( const hhg::Point3D& ) { return 2; };

      x->interpolate( ones, level );
      rhs->interpolate( ones, level );

      hhg::communication::syncFunctionBetweenPrimitives( ( *x ), level );

      P2::face::smoothGSvertexDoF( level,
                                   *face,
                                   p2operator.getVertexToVertexOpr().getFaceStencilID(),
                                   x->getVertexDoFFunction()->getFaceDataID(),
                                   p2operator.getEdgeToVertexOpr().getFaceStencilID(),
                                   x->getEdgeDoFFunction()->getFaceDataID(),
                                   rhs->getVertexDoFFunction()->getFaceDataID() );

      P2::face::smoothGSedgeDoF( level,
                                 *face,
                                 p2operator.getVertexToEdgeOpr().getFaceStencilID(),
                                 x->getVertexDoFFunction()->getFaceDataID(),
                                 p2operator.getEdgeToEdgeOpr().getFaceStencilID(),
                                 x->getEdgeDoFFunction()->getFaceDataID(),
                                 rhs->getEdgeDoFFunction()->getFaceDataID() );

      real_t* edgeDoFData   = face->getData( x->getEdgeDoFFunction()->getFaceDataID() )->getPointer( level );
      real_t* vertexDoFData = face->getData( x->getVertexDoFFunction()->getFaceDataID() )->getPointer( level );

      for( const auto& it : hhg::vertexdof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             vertexDoFData[hhg::vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )],
             1.,
             it.col() << " " << it.row() );
      }

      for( const auto& it : hhg::edgedof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hhg::edgedof::macroface::indexFromVertex( level, it.col(), it.row(), sD::EDGE_HO_E )],
             1.,
             it.col() << " " << it.row() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hhg::edgedof::macroface::indexFromVertex( level, it.col(), it.row(), sD::EDGE_DI_NE )],
             1.,
             it.col() << " " << it.row() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hhg::edgedof::macroface::indexFromVertex( level, it.col(), it.row(), sD::EDGE_VE_N )],
             1.,
             it.col() << " " << it.row() );
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

   P2::edge::smoothGSvertexDoF( level,
                                *doubleEdge,
                                p2operator.getVertexToVertexOpr().getEdgeStencilID(),
                                x->getVertexDoFFunction()->getEdgeDataID(),
                                p2operator.getEdgeToVertexOpr().getEdgeStencilID(),
                                x->getEdgeDoFFunction()->getEdgeDataID(),
                                rhs->getVertexDoFFunction()->getEdgeDataID() );

   P2::edge::smoothGSedgeDoF( level,
                              *doubleEdge,
                              p2operator.getVertexToEdgeOpr().getEdgeStencilID(),
                              x->getVertexDoFFunction()->getEdgeDataID(),
                              p2operator.getEdgeToEdgeOpr().getEdgeStencilID(),
                              x->getEdgeDoFFunction()->getEdgeDataID(),
                              rhs->getEdgeDoFFunction()->getEdgeDataID() );

   real_t* edgeDoFData   = doubleEdge->getData( x->getEdgeDoFFunction()->getEdgeDataID() )->getPointer( level );
   real_t* vertexDoFData = doubleEdge->getData( x->getVertexDoFFunction()->getEdgeDataID() )->getPointer( level );

   for( const auto& it : hhg::vertexdof::macroedge::Iterator( level, 0 ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( vertexDoFData[hhg::vertexdof::macroedge::indexFromVertex( level, it.col(), sD::VERTEX_C )],
                                  1.,
                                  it.col() << " " << it.row() );
   }

   for( const auto& it : hhg::edgedof::macroedge::Iterator( level, 0 ) )
   {
      WALBERLA_CHECK_FLOAT_EQUAL( edgeDoFData[hhg::edgedof::macroedge::indexFromVertex( level, it.col(), sD::EDGE_HO_E )],
                                  1.,
                                  it.col() << " " << it.row() );
   }
}

static void testP2JacobiSmooth()
{
   const uint_t level = 3;

   MeshInfo mesh = MeshInfo::fromGmshFile( "../../data/meshes/quad_2el.msh" );

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

   std::function< real_t( const hhg::Point3D& ) >                          ones = []( const hhg::Point3D& ) { return 13; };
   std::function< real_t( const hhg::Point3D& ) >                          two  = []( const hhg::Point3D& ) { return 2; };
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > onesExtended =
       [&ones]( const hhg::Point3D& xx, const std::vector< real_t >& ) { return ones( xx ); };

   tmp->interpolate( ones, level );
   rhs->interpolate( ones, level );

   hhg::communication::syncFunctionBetweenPrimitives( ( *tmp ), level );
   hhg::communication::syncFunctionBetweenPrimitives( ( *rhs ), level );

   for( auto e : storage->getEdges() )
   {
      Edge* edge = e.second.get();
      vertexdof::macroedge::interpolate( level, *edge, x->getVertexDoFFunction()->getEdgeDataID(), {}, onesExtended );
      edgedof::macroedge::interpolate( level, *edge, x->getEdgeDoFFunction()->getEdgeDataID(), {}, onesExtended );
   }
   for( auto v : storage->getVertices() )
   {
      Vertex* vertex = v.second.get();
      vertexdof::macrovertex::interpolate( *vertex, x->getVertexDoFFunction()->getVertexDataID(), {}, onesExtended, level );
   }

   hhg::communication::syncFunctionBetweenPrimitives( ( *x ), level );

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
                                            tmp->getVertexDoFFunction()->getFaceDataID(),
                                            x->getVertexDoFFunction()->getFaceDataID(),
                                            p2operator.getEdgeToVertexOpr().getFaceStencilID(),
                                            tmp->getEdgeDoFFunction()->getFaceDataID(),
                                            rhs->getVertexDoFFunction()->getFaceDataID() );

      P2::macroface::smoothJacobiEdgeDoF( level,
                                          *face,
                                          p2operator.getVertexToEdgeOpr().getFaceStencilID(),
                                          tmp->getVertexDoFFunction()->getFaceDataID(),
                                          p2operator.getEdgeToEdgeOpr().getFaceStencilID(),
                                          tmp->getEdgeDoFFunction()->getFaceDataID(),
                                          x->getEdgeDoFFunction()->getFaceDataID(),
                                          rhs->getEdgeDoFFunction()->getFaceDataID() );

      real_t* edgeDoFData   = face->getData( x->getEdgeDoFFunction()->getFaceDataID() )->getPointer( level );
      real_t* vertexDoFData = face->getData( x->getVertexDoFFunction()->getFaceDataID() )->getPointer( level );

      for( const auto& it : hhg::vertexdof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             vertexDoFData[hhg::vertexdof::macroface::indexFromVertex( level, it.col(), it.row(), sD::VERTEX_C )],
             13.,
             it.col() << " " << it.row() );
      }

      for( const auto& it : hhg::edgedof::macroface::Iterator( level, 0 ) )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hhg::edgedof::macroface::indexFromVertex( level, it.col(), it.row(), sD::EDGE_HO_E )],
             13.,
             it.col() << " " << it.row() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hhg::edgedof::macroface::indexFromVertex( level, it.col(), it.row(), sD::EDGE_DI_NE )],
             13.,
             it.col() << " " << it.row() );

         WALBERLA_CHECK_FLOAT_EQUAL(
             edgeDoFData[hhg::edgedof::macroface::indexFromVertex( level, it.col(), it.row(), sD::EDGE_VE_N )],
             13.,
             it.col() << " " << it.row() );
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

   P2::edge::smoothGSvertexDoF( level,
                                *doubleEdge,
                                p2operator.getVertexToVertexOpr().getEdgeStencilID(),
                                x->getVertexDoFFunction()->getEdgeDataID(),
                                p2operator.getEdgeToVertexOpr().getEdgeStencilID(),
                                x->getEdgeDoFFunction()->getEdgeDataID(),
                                rhs->getVertexDoFFunction()->getEdgeDataID() );

   P2::edge::smoothGSedgeDoF( level,
                              *doubleEdge,
                              p2operator.getVertexToEdgeOpr().getEdgeStencilID(),
                              x->getVertexDoFFunction()->getEdgeDataID(),
                              p2operator.getEdgeToEdgeOpr().getEdgeStencilID(),
                              x->getEdgeDoFFunction()->getEdgeDataID(),
                              rhs->getEdgeDoFFunction()->getEdgeDataID() );

   ///TODO: enable once jacobi on macroedges is implemented
#if 0
  real_t *edgeDoFData = doubleEdge->getData(x->getEdgeDoFFunction()->getEdgeDataID())->getPointer(level);
  real_t *vertexDoFData = doubleEdge->getData(x->getVertexDoFFunction()->getEdgeDataID())->getPointer(level);

  for (const auto &it : hhg::vertexdof::macroedge::Iterator(level, 0)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      vertexDoFData[hhg::vertexdof::macroedge::indexFromVertex(level,it.col(), sD::VERTEX_C)],
      1.,
      it.col() << " " << it.row());
  }

  for (const auto &it : hhg::edgedof::macroedge::Iterator(level, 0)) {
    WALBERLA_CHECK_FLOAT_EQUAL(
      edgeDoFData[hhg::edgedof::macroedge::indexFromVertex(level,it.col(), sD::EDGE_HO_E)],
      1.,
      it.col() << " " << it.row());
  }
#endif
}

} // namespace hhg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testP2Smooth();
   hhg::testP2JacobiSmooth();

   return EXIT_SUCCESS;
}