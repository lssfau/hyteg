#include "P2ConstantOperator.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-conversion"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "generated/p2_diffusion.h"
#include "generated/p2_div.h"
#include "generated/p2_divt.h"
#include "generated/p2_mass.h"
#include "generated/p2_tet_diffusion.h"

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "tinyhhg_core/p2functionspace/P2Elements.hpp"
#include "tinyhhg_core/p2functionspace/P2MacroEdge.hpp"
#include "tinyhhg_core/p2functionspace/P2MacroFace.hpp"
#include "tinyhhg_core/p2functionspace/P2Smooth.hpp"

namespace hhg {

template < class UFCOperator2D, class UFCOperator3D >
P2ConstantOperator< UFCOperator2D, UFCOperator3D >::P2ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                        size_t                                     minLevel,
                                                                        size_t                                     maxLevel )
: Operator( storage, minLevel, maxLevel )
, vertexToVertex( storage, minLevel, maxLevel )
, edgeToVertex( storage, minLevel, maxLevel )
, vertexToEdge( storage, minLevel, maxLevel )
, edgeToEdge( storage, minLevel, maxLevel )
{
   if( storage_->hasGlobalCells() )
   {
      const bool assemblyDefined = !std::is_same< UFCOperator3D, hhg::fenics::UndefinedAssembly >::value;
      WALBERLA_CHECK( assemblyDefined, "Assembly undefined for 3D elements." );
      if( !std::is_same< UFCOperator3D, fenics::NoAssemble >::value )
      {
         assembleStencils3D();
      }
   } else
   {
      if( !std::is_same< UFCOperator2D, fenics::NoAssemble >::value )
      {
         const bool assemblyDefined = !std::is_same< UFCOperator2D, hhg::fenics::UndefinedAssembly >::value;
         WALBERLA_CHECK( assemblyDefined, "Assembly undefined for 2D elements." );
         assembleStencils();
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::assembleStencils()
{
   using namespace P2Elements;

   // Initialize memory for local 6x6 matrices
   Matrix6r local_stiffness_gray;
   Matrix6r local_stiffness_blue;

   // Assemble stencils on all levels
   for( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      // Assemble face stencils
      for( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         // Compute both local stiffness matrices
         compute_local_stiffness( face, level, local_stiffness_gray, fenics::GRAY );
         compute_local_stiffness( face, level, local_stiffness_blue, fenics::BLUE );

         //        WALBERLA_LOG_DEVEL_ON_ROOT("local_stiffness_gray =\n" << local_stiffness_gray);
         //        WALBERLA_LOG_DEVEL_ON_ROOT("local_stiffness_blue =\n" << local_stiffness_blue);

         // Assemble vertexToVertex stencil
         real_t* vStencil = storage_->getFace( face.getID() )->getData( vertexToVertex.getFaceStencilID() )->getPointer( level );
         P2Face::VertexToVertex::assembleStencil( local_stiffness_gray, local_stiffness_blue, vStencil );
         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Face = {}", PointND<real_t, 7>(&vStencil[0])));

         // Assemble edgeToVertex stencil
         vStencil = storage_->getFace( face.getID() )->getData( edgeToVertex.getFaceStencilID() )->getPointer( level );
         P2Face::EdgeToVertex::assembleStencil( local_stiffness_gray, local_stiffness_blue, vStencil );
         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Face = {}", PointND<real_t, 12>(&vStencil[0])));

         // Assemble vertexToEdge stencil
         vStencil = storage_->getFace( face.getID() )->getData( vertexToEdge.getFaceStencilID() )->getPointer( level );
         P2Face::VertexToEdge::assembleStencil( local_stiffness_gray, local_stiffness_blue, vStencil );
         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Face = {}", PointND<real_t, 12>(&vStencil[0])));

         // Assemble edgeToEdge stencil
         vStencil = storage_->getFace( face.getID() )->getData( edgeToEdge.getFaceStencilID() )->getPointer( level );
         P2Face::EdgeToEdge::assembleStencil( local_stiffness_gray, local_stiffness_blue, vStencil );
         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToEdge/Face = {}", PointND<real_t, 15>(&vStencil[0])));
      }

      // Assemble edge stencils
      for( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         // Assemble vertexToVertex stencil
         Face*   face     = storage_->getFace( edge.neighborFaces()[0] );
         real_t* vStencil = storage_->getEdge( edge.getID() )->getData( vertexToVertex.getEdgeStencilID() )->getPointer( level );
         compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
         compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
         P2Edge::VertexToVertex::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true );

         if( edge.getNumNeighborFaces() == 2 )
         {
            face = storage_->getFace( edge.neighborFaces()[1] );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
            P2Edge::VertexToVertex::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false );
         }

         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Edge = {}", PointND<real_t, 7>(&vStencil[0])));

         // Assemble edgeToVertex
         face     = storage_->getFace( edge.neighborFaces()[0] );
         vStencil = storage_->getEdge( edge.getID() )->getData( edgeToVertex.getEdgeStencilID() )->getPointer( level );
         compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
         compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
         P2Edge::EdgeToVertex::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true );

         if( edge.getNumNeighborFaces() == 2 )
         {
            face = storage_->getFace( edge.neighborFaces()[1] );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
            P2Edge::EdgeToVertex::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false );
         }

         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Edge = {}", PointND<real_t, 7>(&vStencil[0])));

         // Assemble vertexToEdge stencil
         face     = storage_->getFace( edge.neighborFaces()[0] );
         vStencil = storage_->getEdge( edge.getID() )->getData( vertexToEdge.getEdgeStencilID() )->getPointer( level );
         compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
         compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
         P2Edge::VertexToEdge::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true );

         if( edge.getNumNeighborFaces() == 2 )
         {
            face = storage_->getFace( edge.neighborFaces()[1] );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
            P2Edge::VertexToEdge::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false );
         }

         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Edge = {}", PointND<real_t, 4>(&vStencil[0])));

         // Assemble edgeToEdge stencil
         face     = storage_->getFace( edge.neighborFaces()[0] );
         vStencil = storage_->getEdge( edge.getID() )->getData( edgeToEdge.getEdgeStencilID() )->getPointer( level );
         compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
         compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
         P2Edge::EdgeToEdge::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true );

         if( edge.getNumNeighborFaces() == 2 )
         {
            face = storage_->getFace( edge.neighborFaces()[1] );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
            P2Edge::EdgeToEdge::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false );
         }

         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToEdge/Edge = {}", PointND<real_t, 5>(&vStencil[0])));
      }

      for( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         // Assemble VertexToVertex
         real_t* vStencil =
             storage_->getVertex( vertex.getID() )->getData( vertexToVertex.getVertexStencilID() )->getPointer( level );
         for( auto& faceId : vertex.neighborFaces() )
         {
            Face* face = storage_->getFace( faceId );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            P2Vertex::VertexToVertex::assembleStencil( vertex, *face, local_stiffness_gray, vStencil, storage_ );
         }

         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Vertex = {}", PointND<real_t, 5>(&vStencil[0])));

         // Assemble EdgeToVertex
         vStencil = storage_->getVertex( vertex.getID() )->getData( edgeToVertex.getVertexStencilID() )->getPointer( level );
         for( auto& faceId : vertex.neighborFaces() )
         {
            Face* face = storage_->getFace( faceId );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            P2Vertex::EdgeToVertex::assembleStencil( vertex, *face, local_stiffness_gray, vStencil, storage_ );
         }

         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Vertex = {}", PointND<real_t, 5>(&vStencil[0])));
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::assembleStencils3D()
{
   UFCOperator3D ufcOperator;

   // Assemble stencils on all levels
   for( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      /////////////
      /// Faces ///
      /////////////

      for( const auto& it : storage_->getFaces() )
      {
         const auto& face = *it.second;

         WALBERLA_ASSERT_GREATER( face.getNumNeighborCells(), 0 );

         ////////////////////////
         /// Vertex -> Vertex ///
         ////////////////////////

         auto stencilMemory = face.getData( getVertexToVertexOpr().getFaceStencilID() )->getPointer( level );

         auto stencil =
             P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, face, indexing::Index( 1, 1, 0 ), level, ufcOperator );

         if( face.getNumNeighborCells() == 1 )
         {
            for( const auto stencilDir : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter )
            {
               if( stencil.count( stencilDir ) == 0 )
               {
                  stencil[stencilDir] = real_c( 0 );
               }
            }
         } else
         {
            for( const auto stencilDir : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter )
            {
               if( stencil.count( stencilDir ) == 0 )
               {
                  stencil[stencilDir] = real_c( 0 );
               }
            }
         }

         for( const auto stencilIt : stencil )
         {
            const auto stencilIdx     = vertexdof::stencilIndexFromVertex( stencilIt.first );
            stencilMemory[stencilIdx] = stencil[stencilIt.first];
         }

         for( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
         {
            const auto& cell = *( storage_->getCell( face.neighborCells().at( neighborCellID ) ) );

            const std::array< edgedof::EdgeDoFOrientation, 3 > faceOrientations = {
                edgedof::EdgeDoFOrientation::X, edgedof::EdgeDoFOrientation::Y, edgedof::EdgeDoFOrientation::XY};

            const uint_t                  localFaceID          = cell.getLocalFaceID( face.getID() );
            const std::array< uint_t, 4 > localVertexIDsAtCell = {
                cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
                cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
                cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ),
                6 - cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ) -
                    cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ) -
                    cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 )};
            const auto vertexAssemblyIndexInCell = indexing::basisConversion(
                indexing::Index( 1, 1, 0 ), localVertexIDsAtCell, {0, 1, 2, 3}, levelinfo::num_microvertices_per_edge( level ) );
            const auto edgeAssemblyIndexInCell = indexing::basisConversion(
                indexing::Index( 1, 1, 0 ), localVertexIDsAtCell, {0, 1, 2, 3}, levelinfo::num_microedges_per_edge( level ) );

            //////////////////////
            /// Edge -> Vertex ///
            //////////////////////

            auto& edgeToVertexStencilMemory = face.getData( getEdgeToVertexOpr().getFaceStencil3DID() )->getData( level );
            for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
               edgeToVertexStencilMemory[neighborCellID][leafOrientation] =
                   P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                       vertexAssemblyIndexInCell, leafOrientation, cell, level, ufcOperator );
            }

            //////////////////////
            /// Vertex -> Edge ///
            //////////////////////

            auto& vertexToEdgeStencilMemory = face.getData( getVertexToEdgeOpr().getFaceStencil3DID() )->getData( level );
            for( const auto& centerOrientation : faceOrientations )
            {
               const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientation(
                   centerOrientation, localVertexIDsAtCell.at( 0 ), localVertexIDsAtCell.at( 1 ), localVertexIDsAtCell.at( 2 ) );
               vertexToEdgeStencilMemory[neighborCellID][convertedCenterOrientation] =
                   P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                       edgeAssemblyIndexInCell, convertedCenterOrientation, cell, level, ufcOperator );
            }

            ////////////////////
            /// Edge -> Edge ///
            ////////////////////

            auto& edgeToEdgeStencilMemory = face.getData( getEdgeToEdgeOpr().getFaceStencil3DID() )->getData( level );
            for( const auto& centerOrientation : faceOrientations )
            {
               for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
               {
                  const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientation( centerOrientation,
                                                                                              localVertexIDsAtCell.at( 0 ),
                                                                                              localVertexIDsAtCell.at( 1 ),
                                                                                              localVertexIDsAtCell.at( 2 ) );
                  edgeToEdgeStencilMemory[neighborCellID][convertedCenterOrientation][leafOrientation] =
                      P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                          edgeAssemblyIndexInCell, convertedCenterOrientation, leafOrientation, cell, level, ufcOperator );
               }
            }
         }
      }

      /////////////
      /// Cells ///
      /////////////

      for( const auto& it : storage_->getCells() )
      {
         const auto& cell = *it.second;

         // vertex to vertex
         auto       vertexToVertexStencilMemory = cell.getData( getVertexToVertexOpr().getCellStencilID() )->getPointer( level );
         const auto vertexToVertexStencilMap    = P1Elements::P1Elements3D::assembleP1LocalStencil(
             getStorage(), cell, indexing::Index( 1, 1, 1 ), level, ufcOperator );

         for( const auto stencilIt : vertexToVertexStencilMap )
         {
            const auto stencilIdx                   = vertexdof::stencilIndexFromVertex( stencilIt.first );
            vertexToVertexStencilMemory[stencilIdx] = stencilIt.second;
         }

         // edge to vertex
         auto& edgeToVertexStencilMemory = cell.getData( getEdgeToVertexOpr().getCellStencilID() )->getData( level );
         for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
         {
            const auto edgeToVertexStencilMap = P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                indexing::Index( 1, 1, 1 ), leafOrientation, cell, level, ufcOperator );
            for( const auto stencilIt : edgeToVertexStencilMap )
            {
               edgeToVertexStencilMemory[leafOrientation][stencilIt.first] = stencilIt.second;
            }
         }

         // vertex to edge
         auto& vertexToEdgeStencilMemory = cell.getData( getVertexToEdgeOpr().getCellStencilID() )->getData( level );
         for( const auto& centerOrientation : edgedof::allEdgeDoFOrientations )
         {
            const auto vertexToEdgeStencilMap = P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                edgedof::macrocell::getInnerIndexByOrientation( centerOrientation ),
                centerOrientation,
                cell,
                level,
                ufcOperator );
            for( const auto stencilIt : vertexToEdgeStencilMap )
            {
               vertexToEdgeStencilMemory[centerOrientation][stencilIt.first] = stencilIt.second;
            }
         }

         // edge to edge
         auto& edgeToEdgeStencilMemory = cell.getData( getEdgeToEdgeOpr().getCellStencilID() )->getData( level );
         for( const auto& centerOrientation : edgedof::allEdgeDoFOrientations )
         {
            for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
               const auto edgeToEdgeStencilMap = P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                   edgedof::macrocell::getInnerIndexByOrientation( centerOrientation ),
                   centerOrientation,
                   leafOrientation,
                   cell,
                   level,
                   ufcOperator );
               for( const auto stencilIt : edgeToEdgeStencilMap )
               {
                  edgeToEdgeStencilMemory[centerOrientation][leafOrientation][stencilIt.first] = stencilIt.second;
               }
            }
         }
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::apply_impl( P2Function< real_t >& src,
                                                                     P2Function< real_t >& dst,
                                                                     size_t                level,
                                                                     DoFType               flag,
                                                                     UpdateType            updateType )
{
   vertexToVertex.apply( *src.getVertexDoFFunction(), *dst.getVertexDoFFunction(), level, flag, updateType );
   edgeToVertex.apply( *src.getEdgeDoFFunction(), *dst.getVertexDoFFunction(), level, flag, Add );

   edgeToEdge.apply( *src.getEdgeDoFFunction(), *dst.getEdgeDoFFunction(), level, flag, updateType );
   vertexToEdge.apply( *src.getVertexDoFFunction(), *dst.getEdgeDoFFunction(), level, flag, Add );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::smooth_gs_impl( P2Function< real_t >& dst,
                                                                         P2Function< real_t >& rhs,
                                                                         size_t                level,
                                                                         DoFType               flag )
{
   dst.getVertexDoFFunction()->communicate< Face, Edge >( level );
   dst.getVertexDoFFunction()->communicate< Edge, Vertex >( level );
   dst.getEdgeDoFFunction()->communicate< Face, Edge >( level );
   dst.getEdgeDoFFunction()->communicate< Edge, Vertex >( level );

   for( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         P2::vertex::smoothGSvertexDoF( level,
                                        vertex,
                                        vertexToVertex.getVertexStencilID(),
                                        dst.getVertexDoFFunction()->getVertexDataID(),
                                        edgeToVertex.getVertexStencilID(),
                                        dst.getEdgeDoFFunction()->getVertexDataID(),
                                        rhs.getVertexDoFFunction()->getVertexDataID() );
      }
   }

   dst.getVertexDoFFunction()->communicate< Vertex, Edge >( level );
   dst.getEdgeDoFFunction()->communicate< Vertex, Edge >( level );

   for( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::smoothGaussSeidel( level,
                                           edge,
                                           vertexToVertex.getEdgeStencilID(),
                                           edgeToVertex.getEdgeStencilID(),
                                           dst.getVertexDoFFunction()->getEdgeDataID(),
                                           vertexToEdge.getEdgeStencilID(),
                                           edgeToEdge.getEdgeStencilID(),
                                           dst.getEdgeDoFFunction()->getEdgeDataID(),
                                           rhs.getVertexDoFFunction()->getEdgeDataID(),
                                           rhs.getEdgeDoFFunction()->getEdgeDataID() );
      }
   }

   dst.getVertexDoFFunction()->communicate< Edge, Face >( level );
   dst.getEdgeDoFFunction()->communicate< Edge, Face >( level );

   for( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         P2::macroface::smoothGaussSeidel( level,
                                           face,
                                           vertexToVertex.getFaceStencilID(),
                                           edgeToVertex.getFaceStencilID(),
                                           dst.getVertexDoFFunction()->getFaceDataID(),
                                           vertexToEdge.getFaceStencilID(),
                                           edgeToEdge.getFaceStencilID(),
                                           dst.getEdgeDoFFunction()->getFaceDataID(),
                                           rhs.getVertexDoFFunction()->getFaceDataID(),
                                           rhs.getEdgeDoFFunction()->getFaceDataID() );
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::smooth_jac_impl( P2Function< real_t >& dst,
                                                                          P2Function< real_t >& rhs,
                                                                          P2Function< real_t >& src,
                                                                          size_t                level,
                                                                          DoFType               flag )
{
   ///TODO: remove unneccessary communication here
   src.getVertexDoFFunction()->communicate< Face, Edge >( level );
   src.getVertexDoFFunction()->communicate< Edge, Vertex >( level );
   src.getVertexDoFFunction()->communicate< Vertex, Edge >( level );
   src.getVertexDoFFunction()->communicate< Edge, Face >( level );
   src.getEdgeDoFFunction()->communicate< Face, Edge >( level );
   src.getEdgeDoFFunction()->communicate< Edge, Vertex >( level );
   src.getEdgeDoFFunction()->communicate< Vertex, Edge >( level );
   src.getEdgeDoFFunction()->communicate< Edge, Face >( level );

   for( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         P2::macroface::smoothJacobiVertexDoF( level,
                                               face,
                                               vertexToVertex.getFaceStencilID(),
                                               src.getVertexDoFFunction()->getFaceDataID(),
                                               dst.getVertexDoFFunction()->getFaceDataID(),
                                               edgeToVertex.getFaceStencilID(),
                                               src.getEdgeDoFFunction()->getFaceDataID(),
                                               rhs.getVertexDoFFunction()->getFaceDataID() );
      }
   }
   for( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         P2::macroface::smoothJacobiEdgeDoF( level,
                                             face,
                                             vertexToEdge.getFaceStencilID(),
                                             src.getVertexDoFFunction()->getFaceDataID(),
                                             edgeToEdge.getFaceStencilID(),
                                             src.getEdgeDoFFunction()->getFaceDataID(),
                                             dst.getEdgeDoFFunction()->getFaceDataID(),
                                             rhs.getEdgeDoFFunction()->getFaceDataID() );
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::compute_local_stiffness( const Face&         face,
                                                                                  size_t              level,
                                                                                  Matrix6r&           local_stiffness,
                                                                                  fenics::ElementType element_type )
{
   real_t coords[6];
   fenics::compute_micro_coords( face, level, coords, element_type );
   UFCOperator2D gen;
   gen.tabulate_tensor( local_stiffness.data(), NULL, coords, 0 );
}

template class P2ConstantOperator< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise >;
template class P2ConstantOperator< p2_mass_cell_integral_0_otherwise >;

template class P2ConstantOperator< p2_divt_cell_integral_0_otherwise >;
template class P2ConstantOperator< p2_divt_cell_integral_1_otherwise >;
template class P2ConstantOperator< p2_div_cell_integral_0_otherwise >;
template class P2ConstantOperator< p2_div_cell_integral_1_otherwise >;

} // namespace hhg