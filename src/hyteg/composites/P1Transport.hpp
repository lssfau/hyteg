#pragma once

#include "hyteg/composites/transport/VertexDoFMacroCellTransport.hpp"
#include "hyteg/composites/transport/VertexDoFMacroEdgeTransport.hpp"
#include "hyteg/composites/transport/VertexDoFMacroFaceTransport.hpp"
#include "hyteg/composites/transport/VertexDoFMacroVertexTransport.hpp"

namespace hyteg {

class P1Transport
{
 public:
   P1Transport( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t minLevel, const uint_t maxLevel )
   : storage_( storage )
   , invLumpedMass_( storage, minLevel, maxLevel )
   , A_( storage, minLevel, maxLevel )
   , divT_x_( storage, minLevel, maxLevel )
   , divT_y_( storage, minLevel, maxLevel )
   , divT_z_( storage, minLevel, maxLevel )
   , tmp0_( "tmp0", storage, minLevel, maxLevel )
   , tmp1_( "tmp1", storage, minLevel, maxLevel )
   , tmp2_( "tmp2", storage, minLevel, maxLevel )
   , tmp3_( "tmp3", storage, minLevel, maxLevel )
   , tmp4_( "tmp4", storage, minLevel, maxLevel )
   , interSol_( "interSol", storage, minLevel, maxLevel )
   , cOld_( "cOld", storage, minLevel, maxLevel )
   {}

   void step( P1Function< real_t >& c,
              P1Function< real_t >& ux,
              P1Function< real_t >& uy,
              P1Function< real_t >& uz,
              const uint_t&         level,
              const DoFType&        flag,
              const real_t&         dt,
              const real_t&         a )
   {
      cOld_.assign( {1.0}, {c}, level, flag );

      // diffusive term
      A_.apply( c, tmp0_, level, flag, Replace );
      tmp0_.assign( {a * dt}, {tmp0_}, level, flag );

      //      transportStupidApply( c, tmp4_, ux, uy, uz, level, flag );
      transportApply< true >( c, tmp4_, ux, uy, uz, level, flag );

      tmp0_.add( {-dt}, {tmp4_}, level, flag );
      invLumpedMass_.apply( tmp0_, c, level, flag, Add );

      //    interSol_.assign({0.5, 0.5}, {&c, &cOld_}, level, flag);
      //
      //    A_.apply( c, tmp0_, level, flag, Replace );
      //    tmp0_.assign( {a * 0.5 * dt}, {&tmp0_}, level, flag );
      //
      //    divT_x_.apply( c, tmp1_, level, flag, Replace );
      //    divT_y_.apply( c, tmp2_, level, flag, Replace );
      //    divT_z_.apply( c, tmp3_, level, flag, Replace );
      //    tmp1_.multElementwise( {&ux, &tmp1_}, level, flag );
      //    tmp2_.multElementwise( {&uy, &tmp2_}, level, flag );
      //    tmp3_.multElementwise( {&uz, &tmp3_}, level, flag );
      //
      //    tmp4_.assign( {0.5 * dt, 0.5 * dt, 0.5 * dt}, {&tmp1_, &tmp2_, &tmp3_}, level, flag );
      //
      //    tmp0_.add( {-1.0}, {&tmp4_}, level, flag );
      //    invLumpedMass_.apply( tmp0_, interSol_, level, flag, Add );

      //    for (uint_t i = 0; i < 3; ++i)
      //    {
      //      A_.apply(c, tmp0_, level, flag);
      //      A_.smooth_gs(c, tmp0_, level, flag);
      //    }
   }

 private:
   void transportStupidApply( P1Function< real_t >& in,
                              P1Function< real_t >& out,
                              P1Function< real_t >& ux,
                              P1Function< real_t >& uy,
                              P1Function< real_t >& uz,
                              const uint_t&         level,
                              const DoFType&        flag )
   {
      divT_x_.apply( in, tmp1_, level, flag, Replace );
      divT_y_.apply( in, tmp2_, level, flag, Replace );
      divT_z_.apply( in, tmp3_, level, flag, Replace );
      tmp1_.multElementwise( {&ux, &tmp1_}, level, flag );
      tmp2_.multElementwise( {&uy, &tmp2_}, level, flag );
      tmp3_.multElementwise( {&uz, &tmp3_}, level, flag );

      out.assign( {1.0, 1.0, 1.0}, {tmp1_, tmp2_, tmp3_}, level, flag );
   }

   template < bool AlgebraicUpwind >
   void transportApply( P1Function< real_t >& src,
                        P1Function< real_t >& dst,
                        P1Function< real_t >& ux,
                        P1Function< real_t >& uy,
                        P1Function< real_t >& uz,
                        const uint_t&         level,
                        const DoFType&        flag )
   {
      src.communicate< Vertex, Edge >( level );
      src.communicate< Edge, Face >( level );
      src.communicate< Face, Cell >( level );

      src.communicate< Cell, Face >( level );
      src.communicate< Face, Edge >( level );
      src.communicate< Edge, Vertex >( level );

      for( const auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            vertexdof::transport::macrovertex::apply< real_t, AlgebraicUpwind >( level,
                                                                                 vertex,
                                                                                 src.getVertexDataID(),
                                                                                 dst.getVertexDataID(),
                                                                                 ux.getVertexDataID(),
                                                                                 uy.getVertexDataID(),
                                                                                 uz.getVertexDataID(),
                                                                                 divT_x_.getVertexStencilID(),
                                                                                 divT_y_.getVertexStencilID(),
                                                                                 divT_z_.getVertexStencilID() );
         }
      }

      for( const auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            vertexdof::transport::macroedge::apply< real_t, AlgebraicUpwind >( level,
                                                                               edge,
                                                                               src.getEdgeDataID(),
                                                                               dst.getEdgeDataID(),
                                                                               ux.getEdgeDataID(),
                                                                               uy.getEdgeDataID(),
                                                                               uz.getEdgeDataID(),
                                                                               divT_x_.getEdgeStencilID(),
                                                                               divT_y_.getEdgeStencilID(),
                                                                               divT_z_.getEdgeStencilID() );
         }
      }

      for( const auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            vertexdof::transport::macroface::apply< real_t, AlgebraicUpwind >( level,
                                                                               face,
                                                                               src.getFaceDataID(),
                                                                               dst.getFaceDataID(),
                                                                               ux.getFaceDataID(),
                                                                               uy.getFaceDataID(),
                                                                               uz.getFaceDataID(),
                                                                               divT_x_.getFaceStencilID(),
                                                                               divT_y_.getFaceStencilID(),
                                                                               divT_z_.getFaceStencilID() );
         }
      }

      for( const auto& it : storage_->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if( testFlag( cellBC, flag ) )
         {
            vertexdof::transport::macrocell::apply< real_t, AlgebraicUpwind >( level,
                                                                               cell,
                                                                               src.getCellDataID(),
                                                                               dst.getCellDataID(),
                                                                               ux.getCellDataID(),
                                                                               uy.getCellDataID(),
                                                                               uz.getCellDataID(),
                                                                               divT_x_.getCellStencilID(),
                                                                               divT_y_.getCellStencilID(),
                                                                               divT_z_.getCellStencilID() );
         }
      }
   }

   const std::shared_ptr< PrimitiveStorage > storage_;
   P1LumpedInvMassOperator                   invLumpedMass_;
   P1ConstantLaplaceOperator                 A_;
   P1DivTxOperator                           divT_x_;
   P1DivTyOperator                           divT_y_;
   P1DivTzOperator                           divT_z_;
   P1Function< real_t >                      tmp0_;
   P1Function< real_t >                      tmp1_;
   P1Function< real_t >                      tmp2_;
   P1Function< real_t >                      tmp3_;
   P1Function< real_t >                      tmp4_;
   P1Function< real_t >                      interSol_;
   P1Function< real_t >                      cOld_;
};

} // namespace hyteg
