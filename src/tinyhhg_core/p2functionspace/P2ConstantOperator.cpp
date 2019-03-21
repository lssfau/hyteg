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
#include "generated/p2_tet_diffusion.h"
#include "generated/p2_tet_mass.h"

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

template < class P2Form >
P2ConstantOperator< P2Form >::P2ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                        size_t                                     minLevel,
                                                                        size_t                                     maxLevel )
: Operator( storage, minLevel, maxLevel )
, vertexToVertex( storage, minLevel, maxLevel )
, edgeToVertex( storage, minLevel, maxLevel )
, vertexToEdge( storage, minLevel, maxLevel )
, edgeToEdge( storage, minLevel, maxLevel )
{
}

template < class P2Form >
void P2ConstantOperator< P2Form >::apply(const P2Function< real_t >& src,
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

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_gs(const P2Function <real_t> &dst,
                                                                   const P2Function <real_t> &rhs,
                                                                   const size_t level,
                                                                   const DoFType flag) const
{
  smooth_sor( dst, rhs, 1.0, level, flag );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_sor(const P2Function <real_t> &dst,
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

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_jac( const P2Function< real_t >& dst,
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

template class P2ConstantOperator< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;

template class P2ConstantOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise > >;

} // namespace hhg
