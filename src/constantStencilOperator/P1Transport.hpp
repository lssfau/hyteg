/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#pragma once

#include "hyteg/composites/transport/VertexDoFMacroCellTransport.hpp"
#include "hyteg/composites/transport/VertexDoFMacroEdgeTransport.hpp"
#include "hyteg/composites/transport/VertexDoFMacroFaceTransport.hpp"
#include "hyteg/composites/transport/VertexDoFMacroVertexTransport.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"

#include "constantStencilOperator/P1ConstantOperator.hpp"

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

   void step( P1Function< real_t >&       c,
              P1VectorFunction< real_t >& vel,
              const uint_t&               level,
              const DoFType&              flag,
              const real_t&               dt,
              const real_t&               a )
   {
      cOld_.assign( { 1.0 }, { c }, level, flag );

      // diffusive term
      A_.apply( c, tmp0_, level, flag, Replace );
      tmp0_.assign( { a * dt }, { tmp0_ }, level, flag );

      //      transportStupidApply( c, tmp4_, ux, uy, uz, level, flag );
      transportApply< true >( c, tmp4_, vel, level, flag );

      tmp0_.add( { -dt }, { tmp4_ }, level, flag );
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
   void transportStupidApply( P1Function< real_t >&       in,
                              P1Function< real_t >&       out,
                              P1VectorFunction< real_t >& vel,
                              const uint_t&               level,
                              const DoFType&              flag )
   {
      divT_x_.apply( in, tmp1_, level, flag, Replace );
      divT_y_.apply( in, tmp2_, level, flag, Replace );
      divT_z_.apply( in, tmp3_, level, flag, Replace );
      tmp1_.multElementwise( { vel[0], tmp1_ }, level, flag );
      tmp2_.multElementwise( { vel[1], tmp2_ }, level, flag );
      if ( vel.getDimension() == 3 )
         tmp3_.multElementwise( { vel[2], tmp3_ }, level, flag );

      out.assign( { 1.0, 1.0, 1.0 }, { tmp1_, tmp2_, tmp3_ }, level, flag );
   }

   template < bool AlgebraicUpwind >
   void transportApply( P1Function< real_t >&       src,
                        P1Function< real_t >&       dst,
                        P1VectorFunction< real_t >& vel,
                        const uint_t&               level,
                        const DoFType&              flag )
   {
      src.communicate< Vertex, Edge >( level );
      src.communicate< Edge, Face >( level );
      src.communicate< Face, Cell >( level );

      src.communicate< Cell, Face >( level );
      src.communicate< Face, Edge >( level );
      src.communicate< Edge, Vertex >( level );

      for ( const auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            auto uxID = vel[0].getVertexDataID();
            auto uyID = vel[1].getVertexDataID();
            // A hack that is not nice :(
            auto uzID = vel.getDimension() == 3 ? vel[2].getVertexDataID() : vel[0].getVertexDataID();

            vertexdof::transport::macrovertex::apply< real_t, AlgebraicUpwind >( level,
                                                                                 vertex,
                                                                                 src.getVertexDataID(),
                                                                                 dst.getVertexDataID(),
                                                                                 uxID,
                                                                                 uyID,
                                                                                 uzID,
                                                                                 divT_x_.getVertexStencilID(),
                                                                                 divT_y_.getVertexStencilID(),
                                                                                 divT_z_.getVertexStencilID() );
         }
      }

      for ( const auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            auto uxID = vel[0].getEdgeDataID();
            auto uyID = vel[1].getEdgeDataID();
            // A hack that is not nice :(
            auto uzID = vel.getDimension() == 3 ? vel[2].getEdgeDataID() : vel[0].getEdgeDataID();

            vertexdof::transport::macroedge::apply< real_t, AlgebraicUpwind >( level,
                                                                               edge,
                                                                               src.getEdgeDataID(),
                                                                               dst.getEdgeDataID(),
                                                                               uxID,
                                                                               uyID,
                                                                               uzID,
                                                                               divT_x_.getEdgeStencilID(),
                                                                               divT_y_.getEdgeStencilID(),
                                                                               divT_z_.getEdgeStencilID() );
         }
      }

      for ( const auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            auto uxID = vel[0].getFaceDataID();
            auto uyID = vel[1].getFaceDataID();
            // A hack that is not nice :(
            auto uzID = vel.getDimension() == 3 ? vel[2].getFaceDataID() : vel[0].getFaceDataID();

            vertexdof::transport::macroface::apply< real_t, AlgebraicUpwind >( level,
                                                                               face,
                                                                               *src.getStorage(),
                                                                               src.getFaceDataID(),
                                                                               dst.getFaceDataID(),
                                                                               uxID,
                                                                               uyID,
                                                                               uzID,
                                                                               divT_x_.getFaceStencil3DID(),
                                                                               divT_y_.getFaceStencil3DID(),
                                                                               divT_z_.getFaceStencil3DID() );
         }
      }

      for ( const auto& it : storage_->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            auto uxID = vel[0].getCellDataID();
            auto uyID = vel[1].getCellDataID();
            // A hack that is not nice :(
            auto uzID = vel.getDimension() == 3 ? vel[2].getCellDataID() : vel[0].getCellDataID();

            vertexdof::transport::macrocell::apply< real_t, AlgebraicUpwind >( level,
                                                                               cell,
                                                                               src.getCellDataID(),
                                                                               dst.getCellDataID(),
                                                                               uxID,
                                                                               uyID,
                                                                               uzID,
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
