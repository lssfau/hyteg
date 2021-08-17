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

#include <array>

#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/types/pointnd.hpp"

namespace hyteg {

template < class VelocityBaseType >
class DGUpwindOperator : public Operator< DGFunction< real_t >, DGFunction< real_t > >
{
   typedef std::array< VelocityBaseType , 2 > VelocityType;

 public:
   DGUpwindOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                     const VelocityType&                        velocity,
                     size_t                                     minLevel,
                     size_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , velocity_( velocity )
   {}

   ~DGUpwindOperator() override = default;

   void apply( const DGFunction< real_t >& src,
               const DGFunction< real_t >& dst,
               uint_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      // start pulling edge halos
      src.startCommunication< Face, Edge >( level );

      // end pulling edge halos
      src.endCommunication< Face, Edge >( level );

      // start pulling vertex halos
      src.startCommunication< Edge, Vertex >( level );

      // end pulling vertex halos
      src.endCommunication< Edge, Vertex >( level );

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template startCommunication< Edge, Vertex >( level );
      }

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template startCommunication< Face, Edge >( level );
      }

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template endCommunication< Edge, Vertex >( level );
      }

      for( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            dgfunction::macrovertex::upwind< real_t >( level,
                                        vertex,
                                        storage_,
                                        src.getVertexDataID(),
                                        dst.getVertexDataID(),
                                        std::array< PrimitiveDataID< FunctionMemory< real_t >, Vertex >, 2 >{
                                            {velocity_[0].getVertexDataID(), velocity_[1].getVertexDataID()}},
                                        updateType );
         }
      }

      dst.startCommunication< Vertex, Edge >( level );

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template endCommunication< Face, Edge >( level );
      }

      for( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            dgfunction::macroedge::upwind< real_t >( level,
                                      edge,
                                      storage_,
                                      src.getEdgeDataID(),
                                      dst.getEdgeDataID(),
                                      std::array< PrimitiveDataID< FunctionMemory< real_t >, Edge >, 2 >{
                                          {velocity_[0].getEdgeDataID(), velocity_[1].getEdgeDataID()}},
                                      updateType );
         }
      }

      dst.endCommunication< Vertex, Edge >( level );

      dst.startCommunication< Edge, Face >( level );

      for( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            dgfunction::macroface::upwind< real_t >( level,
                                      face,
                                      storage_,
                                      src.getFaceDataID(),
                                      dst.getFaceDataID(),
                                      std::array< PrimitiveDataID< FunctionMemory< real_t >, Face >, 2 >{
                                          {velocity_[0].getFaceDataID(), velocity_[1].getFaceDataID()}},
                                      updateType );
         }
      }

      dst.endCommunication< Edge, Face >( level );
   }

 private:
   VelocityType velocity_;
};

} // namespace hyteg
