#include "P2ConstantOperator.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-conversion"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

#include "generated/p2_diffusion.h"
#include "generated/p2_div.h"
#include "generated/p2_divt.h"
#include "generated/p2_mass.h"
#include "generated/p2_pspg.h"
#include "generated/p2_tet_diffusion.h"
#include "generated/p2_tet_div_tet.h"
#include "generated/p2_tet_divt_tet.h"
#include "generated/p2_tet_mass.h"
#include "generated/p2_tet_pspg_tet.h"

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "tinyhhg_core/p2functionspace/P2Elements.hpp"
#include "tinyhhg_core/p2functionspace/P2MacroEdge.hpp"
#include "tinyhhg_core/p2functionspace/P2MacroFace.hpp"
#include "tinyhhg_core/p2functionspace/P2MacroVertex.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"
#include "tinyhhg_core/p2functionspace/generatedKernels/GeneratedKernelsP2MacroFace2D.hpp"

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
   assembleEdgeToEdgeStencils< UFCOperator3D >( storage_, minLevel_, maxLevel_,
                                                getEdgeToEdgeOpr().getEdgeStencil3DID(),
                                                getEdgeToEdgeOpr().getFaceStencil3DID(),
                                                getEdgeToEdgeOpr().getCellStencilID() );

   assembleVertexToEdgeStencils< UFCOperator3D >( storage_, minLevel_, maxLevel_,
                                                  getVertexToEdgeOpr().getEdgeStencil3DID(),
                                                  getVertexToEdgeOpr().getFaceStencil3DID(),
                                                  getVertexToEdgeOpr().getCellStencilID() );

   assembleEdgeToVertexStencils< UFCOperator3D >( storage_, minLevel_, maxLevel_,
                                                  getEdgeToVertexOpr().getVertexStencil3DID(),
                                                  getEdgeToVertexOpr().getEdgeStencil3DID(),
                                                  getEdgeToVertexOpr().getFaceStencil3DID(),
                                                  getEdgeToVertexOpr().getCellStencilID() );

   UFCOperator3D ufcOperator;

   // Assemble vertex -> vertex stencils on all levels
   for( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      ////////////////
      /// Vertices ///
      ////////////////

      for( const auto& it : storage_->getVertices() )
      {
         const auto & vertex = *it.second;
         WALBERLA_ASSERT_GREATER( vertex.getNumNeighborCells(), 0 );
         auto          stencilSize   = vertex.getData( getVertexToVertexOpr().getVertexStencilID() )->getSize( level );
         auto          stencilMemory = vertex.getData( getVertexToVertexOpr().getVertexStencilID() )->getPointer( level );

         auto stencil = P1Elements::P1Elements3D::assembleP1LocalStencil(
            storage_, vertex, indexing::Index( 0, 0, 0 ), level, ufcOperator );

         WALBERLA_ASSERT_EQUAL( stencilSize, stencil.size() );
         for( uint_t i = 0; i < stencilSize; i++ )
         {
            stencilMemory[i] = stencil[i];
         }
      }

      /////////////
      /// Edges ///
      /////////////

      for( const auto& it : storage_->getEdges() )
      {
         const auto & edge = *it.second;
         WALBERLA_ASSERT_GREATER( edge.getNumNeighborCells(), 0 );

         auto          stencilSize   = edge.getData( getVertexToVertexOpr().getEdgeStencilID() )->getSize( level );
         auto          stencilMemory = edge.getData( getVertexToVertexOpr().getEdgeStencilID() )->getPointer( level );

         auto stencil =
         P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, edge, indexing::Index( 1, 0, 0 ), level, ufcOperator );

         WALBERLA_ASSERT_EQUAL( stencilSize, stencil.size() );
         for( uint_t i = 0; i < stencilSize; i++ )
         {
            stencilMemory[i] = stencil[i];
         }
      }

      /////////////
      /// Faces ///
      /////////////

      for( const auto& it : storage_->getFaces() )
      {
         const auto& face = *it.second;

         WALBERLA_ASSERT_GREATER( face.getNumNeighborCells(), 0 );

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
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::apply(const P2Function< real_t >& src,
                                                               const P2Function< real_t >& dst,
                                                               size_t                level,
                                                               DoFType               flag,
                                                               UpdateType            updateType ) const
{
   vertexToVertex.apply( src.getVertexDoFFunction(), dst.getVertexDoFFunction(), level, flag, updateType );
   edgeToVertex.apply( src.getEdgeDoFFunction(), dst.getVertexDoFFunction(), level, flag, Add );

   edgeToEdge.apply( src.getEdgeDoFFunction(), dst.getEdgeDoFFunction(), level, flag, updateType );
   vertexToEdge.apply( src.getVertexDoFFunction(), dst.getEdgeDoFFunction(), level, flag, Add );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::smooth_gs(const P2Function <real_t> &dst,
                                                                   const P2Function <real_t> &rhs,
                                                                   const size_t level,
                                                                   const DoFType flag) const
{
  smooth_sor( dst, rhs, 1.0, level, flag );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::smooth_sor(const P2Function <real_t> &dst,
                                                                    const P2Function <real_t> &rhs,
                                                                    const real_t& relax,
                                                                    const size_t level,
                                                                    const DoFType flag) const
{
   this->startTiming( "SOR" );

   communication::syncP2FunctionBetweenPrimitives( dst, level );

   for( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::smoothSORVertexDoF( level,
                                              vertex,
                                              relax,
                                              vertexToVertex.getVertexStencilID(),
                                              dst.getVertexDoFFunction().getVertexDataID(),
                                              edgeToVertex.getVertexStencilID(),
                                              dst.getEdgeDoFFunction().getVertexDataID(),
                                              rhs.getVertexDoFFunction().getVertexDataID() );
      }
   }

   dst.getVertexDoFFunction().communicate< Vertex, Edge >( level );
   dst.getEdgeDoFFunction().communicate< Vertex, Edge >( level );

   for( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::smoothSOR( level,
                                   edge,
                                   relax,      
                                   vertexToVertex.getEdgeStencilID(),
                                   edgeToVertex.getEdgeStencilID(),
                                   dst.getVertexDoFFunction().getEdgeDataID(),
                                   vertexToEdge.getEdgeStencilID(),
                                   edgeToEdge.getEdgeStencilID(),
                                   dst.getEdgeDoFFunction().getEdgeDataID(),
                                   rhs.getVertexDoFFunction().getEdgeDataID(),
                                   rhs.getEdgeDoFFunction().getEdgeDataID() );
      }
   }

   dst.getVertexDoFFunction().communicate< Edge, Face >( level );
   dst.getEdgeDoFFunction().communicate< Edge, Face >( level );

   for ( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         if ( globalDefines::useGeneratedKernels && !storage_->hasGlobalCells() )
         {
            real_t* v_dst_data = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
            real_t* v_rhs_data = face.getData( rhs.getVertexDoFFunction().getFaceDataID() )->getPointer( level );

            real_t* e_dst_data = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
            real_t* e_rhs_data = face.getData( rhs.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

            real_t* v2v_opr_data = face.getData( vertexToVertex.getFaceStencilID() )->getPointer( level );
            real_t* v2e_opr_data = face.getData( vertexToEdge.getFaceStencilID() )->getPointer( level );
            real_t* e2v_opr_data = face.getData( edgeToVertex.getFaceStencilID() )->getPointer( level );
            real_t* e2e_opr_data = face.getData( edgeToEdge.getFaceStencilID() )->getPointer( level );

            P2::macroface::generated::sor_2D_macroface_P2_update_vertexdofs(
                e_dst_data, e2v_opr_data, v_dst_data, v_rhs_data, v2v_opr_data, static_cast< int64_t >( level ), relax );
            P2::macroface::generated::sor_2D_macroface_P2_update_edgedofs( e_dst_data,
                                                                           e_rhs_data,
                                                                           &e2e_opr_data[0],
                                                                           &e2e_opr_data[5],
                                                                           &e2e_opr_data[10],
                                                                           v_dst_data,
                                                                           &v2e_opr_data[0],
                                                                           &v2e_opr_data[4],
                                                                           &v2e_opr_data[8],
                                                                           static_cast< int64_t >( level ),
                                                                           relax );
         }
         else
         {
            P2::macroface::smoothSOR( level,
                                      face,
                                      relax,
                                      vertexToVertex.getFaceStencilID(),
                                      edgeToVertex.getFaceStencilID(),
                                      dst.getVertexDoFFunction().getFaceDataID(),
                                      vertexToEdge.getFaceStencilID(),
                                      edgeToEdge.getFaceStencilID(),
                                      dst.getEdgeDoFFunction().getFaceDataID(),
                                      rhs.getVertexDoFFunction().getFaceDataID(),
                                      rhs.getEdgeDoFFunction().getFaceDataID() );
         }
      }
   }
   this->stopTiming( "SOR" );
}

template < class UFCOperator2D, class UFCOperator3D >
void P2ConstantOperator< UFCOperator2D, UFCOperator3D >::smooth_jac( const P2Function< real_t >& dst,
                                                                     const P2Function< real_t >& rhs,
                                                                     const P2Function< real_t >& src,
                                                                     size_t                      level,
                                                                     DoFType                     flag ) const
{
   ///TODO: remove unneccessary communication here
   src.getVertexDoFFunction().communicate< Face, Edge >( level );
   src.getVertexDoFFunction().communicate< Edge, Vertex >( level );
   src.getVertexDoFFunction().communicate< Vertex, Edge >( level );
   src.getVertexDoFFunction().communicate< Edge, Face >( level );
   src.getEdgeDoFFunction().communicate< Face, Edge >( level );
   src.getEdgeDoFFunction().communicate< Edge, Vertex >( level );
   src.getEdgeDoFFunction().communicate< Vertex, Edge >( level );
   src.getEdgeDoFFunction().communicate< Edge, Face >( level );

   for( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         P2::macroface::smoothJacobiVertexDoF( level,
                                               face,
                                               vertexToVertex.getFaceStencilID(),
                                               src.getVertexDoFFunction().getFaceDataID(),
                                               dst.getVertexDoFFunction().getFaceDataID(),
                                               edgeToVertex.getFaceStencilID(),
                                               src.getEdgeDoFFunction().getFaceDataID(),
                                               rhs.getVertexDoFFunction().getFaceDataID() );
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
                                             src.getVertexDoFFunction().getFaceDataID(),
                                             edgeToEdge.getFaceStencilID(),
                                             src.getEdgeDoFFunction().getFaceDataID(),
                                             dst.getEdgeDoFFunction().getFaceDataID(),
                                             rhs.getEdgeDoFFunction().getFaceDataID() );
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
template class P2ConstantOperator< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise >;

template class P2ConstantOperator< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise >;

template class P2ConstantOperator< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise >;
template class P2ConstantOperator< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise >;
template class P2ConstantOperator< fenics::NoAssemble,                p2_tet_divt_tet_cell_integral_2_otherwise >;
template class P2ConstantOperator< p2_div_cell_integral_0_otherwise,  p2_tet_div_tet_cell_integral_0_otherwise >;
template class P2ConstantOperator< p2_div_cell_integral_1_otherwise,  p2_tet_div_tet_cell_integral_1_otherwise >;
template class P2ConstantOperator< fenics::NoAssemble,                p2_tet_div_tet_cell_integral_2_otherwise >;

} // namespace hhg
