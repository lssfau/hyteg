/*
 * Copyright (c) 2020 Daniel Drzisga.
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
#include "P1ProjectNormalOperator.hpp"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/freeslip/VertexDoFProjectNormal.hpp"

namespace hyteg {

P1ProjectNormalOperator::P1ProjectNormalOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                                                  size_t                                                   minLevel,
                                                  size_t                                                   maxLevel,
                                                  const std::function< void( const Point3D&, Point3D& ) >& normal_function )
: Operator( storage, minLevel, maxLevel )
, normal_function_( normal_function )
{}

void P1ProjectNormalOperator::project( const P1Function< real_t >& dst_u,
                                       const P1Function< real_t >& dst_v,
                                       const P1Function< real_t >& dst_w,
                                       size_t                      level,
                                       DoFType                     flag ) const
{
   this->startTiming( "Project" );
   dst_u.communicate< Vertex, Edge >( level );
   dst_u.communicate< Edge, Face >( level );
   dst_u.communicate< Face, Cell >( level );

   dst_v.communicate< Vertex, Edge >( level );
   dst_v.communicate< Edge, Face >( level );
   dst_v.communicate< Face, Cell >( level );

   if ( storage_->hasGlobalCells() )
   {
      dst_w.communicate< Vertex, Edge >( level );
      dst_w.communicate< Edge, Face >( level );
      dst_w.communicate< Face, Cell >( level );
   }

   dst_u.communicate< Cell, Face >( level );
   dst_u.communicate< Face, Edge >( level );
   dst_u.communicate< Edge, Vertex >( level );

   dst_v.communicate< Cell, Face >( level );
   dst_v.communicate< Face, Edge >( level );
   dst_v.communicate< Edge, Vertex >( level );

   if ( storage_->hasGlobalCells() )
   {
      dst_w.communicate< Cell, Face >( level );
      dst_w.communicate< Face, Edge >( level );
      dst_w.communicate< Edge, Vertex >( level );
   }

   this->timingTree_->start( "Macro-Vertex" );

   for ( const auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst_u.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            vertexdof::macrovertex::projectNormal3D< real_t >( level,
                                                               vertex,
                                                               storage_,
                                                               normal_function_,
                                                               dst_u.getVertexDataID(),
                                                               dst_v.getVertexDataID(),
                                                               dst_w.getVertexDataID() );
         }
         else
         {
            vertexdof::macrovertex::projectNormal2D< real_t >(
                level, vertex, storage_, normal_function_, dst_u.getVertexDataID(), dst_v.getVertexDataID() );
         }
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );

   this->timingTree_->start( "Macro-Edge" );

   if ( level >= 1 )
   {
      for ( const auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst_u.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroedge::projectNormal3D< real_t >(
                   level, edge, storage_, normal_function_, dst_u.getEdgeDataID(), dst_v.getEdgeDataID(), dst_w.getEdgeDataID() );
            }
            else
            {
               vertexdof::macroedge::projectNormal2D< real_t >(
                   level, edge, storage_, normal_function_, dst_u.getEdgeDataID(), dst_v.getEdgeDataID() );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Edge" );

   this->timingTree_->start( "Macro-Face" );

   if ( level >= 2 )
   {
      for ( const auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst_u.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroface::projectNormal3D< real_t >(
                   level, face, storage_, normal_function_, dst_u.getFaceDataID(), dst_v.getFaceDataID(), dst_w.getFaceDataID() );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   this->stopTiming( "Project" );
}

void P1ProjectNormalOperator::project( const P1StokesFunction< real_t >& dst, size_t level, DoFType flag ) const
{
   // This way of delegation will not work for 2D case, as there ist no dst.uvw[2]
   // project( dst.uvw[0], dst.uvw[1], dst.uvw[2], level, flag );

   this->startTiming( "Project" );

   for ( uint_t k = 0; k < dst.uvw().getDimension(); k++ )
   {
      dst.uvw()[k].communicate< Vertex, Edge >( level );
      dst.uvw()[k].communicate< Edge, Face >( level );
      dst.uvw()[k].communicate< Face, Cell >( level );
   }

   for ( uint_t k = 0; k < dst.uvw().getDimension(); k++ )
   {
      dst.uvw()[k].communicate< Cell, Face >( level );
      dst.uvw()[k].communicate< Face, Edge >( level );
      dst.uvw()[k].communicate< Edge, Vertex >( level );
   }

   this->timingTree_->start( "Macro-Vertex" );

   for ( const auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.uvw().getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            vertexdof::macrovertex::projectNormal3D< real_t >( level,
                                                               vertex,
                                                               storage_,
                                                               normal_function_,
                                                               dst.uvw()[0].getVertexDataID(),
                                                               dst.uvw()[1].getVertexDataID(),
                                                               dst.uvw()[2].getVertexDataID() );
         }
         else
         {
            vertexdof::macrovertex::projectNormal2D< real_t >(
                level, vertex, storage_, normal_function_, dst.uvw()[0].getVertexDataID(), dst.uvw()[1].getVertexDataID() );
         }
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );

   this->timingTree_->start( "Macro-Edge" );

   if ( level >= 1 )
   {
      for ( const auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.uvw().getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroedge::projectNormal3D< real_t >( level,
                                                                edge,
                                                                storage_,
                                                                normal_function_,
                                                                dst.uvw()[0].getEdgeDataID(),
                                                                dst.uvw()[1].getEdgeDataID(),
                                                                dst.uvw()[2].getEdgeDataID() );
            }
            else
            {
               vertexdof::macroedge::projectNormal2D< real_t >(
                   level, edge, storage_, normal_function_, dst.uvw()[0].getEdgeDataID(), dst.uvw()[1].getEdgeDataID() );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Edge" );

   this->timingTree_->start( "Macro-Face" );

   if ( level >= 2 )
   {
      for ( const auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.uvw().getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroface::projectNormal3D< real_t >( level,
                                                                face,
                                                                storage_,
                                                                normal_function_,
                                                                dst.uvw()[0].getFaceDataID(),
                                                                dst.uvw()[1].getFaceDataID(),
                                                                dst.uvw()[2].getFaceDataID() );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   this->stopTiming( "Project" );
}

void P1ProjectNormalOperator::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                        const P1Function< idx_t >&                  numU,
                                        const P1Function< idx_t >&                  numV,
                                        const P1Function< idx_t >&                  numW,
                                        uint_t                                      level,
                                        DoFType                                     flag ) const
{
   communication::syncFunctionBetweenPrimitives( numU, level );
   communication::syncFunctionBetweenPrimitives( numV, level );
   communication::syncFunctionBetweenPrimitives( numW, level );

   // The matrix-free application of the projection operator (ID - nn^t) emulates
   // the application of Id by not touching the vector at all.
   // However, the Id diagonal must be assembled.

   for ( const auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = numU.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );

      if ( testFlag( vertexBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            WALBERLA_ABORT( "Sparse matrix assembly not implemented for 3D." );
         }
         else
         {
            vertexdof::macrovertex::saveProjectNormalOperator2D(
                level, vertex, storage_, normal_function_, numU.getVertexDataID(), numV.getVertexDataID(), mat );
         }
      }
      else
      {
         vertexdof::macrovertex::saveIdentityOperator( vertex, numU.getVertexDataID(), mat, level );
         vertexdof::macrovertex::saveIdentityOperator( vertex, numV.getVertexDataID(), mat, level );
         if ( storage_->hasGlobalCells() )
         {
            vertexdof::macrovertex::saveIdentityOperator( vertex, numW.getVertexDataID(), mat, level );
         }
      }
   }

   if ( level >= 1 )
   {
      for ( const auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = numU.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

         if ( testFlag( edgeBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               WALBERLA_ABORT( "Sparse matrix assembly not implemented for 3D." );
            }
            else
            {
               vertexdof::macroedge::saveProjectNormalOperator2D(
                   level, edge, storage_, normal_function_, numU.getEdgeDataID(), numV.getEdgeDataID(), mat );
            }
         }
         else
         {
            vertexdof::macroedge::saveIdentityOperator( level, edge, numU.getEdgeDataID(), mat );
            vertexdof::macroedge::saveIdentityOperator( level, edge, numV.getEdgeDataID(), mat );
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroedge::saveIdentityOperator( level, edge, numW.getEdgeDataID(), mat );
            }
         }
      }
   }

   if ( level >= 2 )
   {
      for ( const auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = numU.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );

         if ( testFlag( faceBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               WALBERLA_ABORT( "Sparse matrix assembly not implemented for 3D." );
            }
            else
            {
               WALBERLA_ABORT( "Normal projection for inner primitives?" );
            }
         }
         else
         {
            vertexdof::macroface::saveIdentityOperator( level, face, numU.getFaceDataID(), mat );
            vertexdof::macroface::saveIdentityOperator( level, face, numV.getFaceDataID(), mat );
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroface::saveIdentityOperator( level, face, numW.getFaceDataID(), mat );
            }
         }
      }
   }
}

} // namespace hyteg
