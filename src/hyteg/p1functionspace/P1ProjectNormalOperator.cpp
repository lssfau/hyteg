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

#include "hyteg/p1functionspace/freeslip/VertexDoFProjectNormal.hpp"

namespace hyteg {

P1ProjectNormalOperator::P1ProjectNormalOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel,
    const std::function<void(const Point3D&, Point3D& )>& normal_function)
: Operator( storage, minLevel, maxLevel ), normal_function_(normal_function)
{
}

void P1ProjectNormalOperator::apply( const P1StokesFunction< real_t >& dst,
                                     size_t                            level,
                                     DoFType                           flag) const
{
   this->startTiming( "Apply" );
   dst.u.communicate< Vertex, Edge >( level );
   dst.u.communicate< Edge, Face >( level );
   dst.u.communicate< Face, Cell >( level );

   dst.v.communicate< Vertex, Edge >( level );
   dst.v.communicate< Edge, Face >( level );
   dst.v.communicate< Face, Cell >( level );

   if ( storage_->hasGlobalCells() )
   {
      dst.w.communicate< Vertex, Edge >( level );
      dst.w.communicate< Edge, Face >( level );
      dst.w.communicate< Face, Cell >( level );
   }

   dst.u.communicate< Cell, Face >( level );
   dst.u.communicate< Face, Edge >( level );
   dst.u.communicate< Edge, Vertex >( level );

   dst.v.communicate< Cell, Face >( level );
   dst.v.communicate< Face, Edge >( level );
   dst.v.communicate< Edge, Vertex >( level );

   if ( storage_->hasGlobalCells() )
   {
      dst.w.communicate< Cell, Face >( level );
      dst.w.communicate< Face, Edge >( level );
      dst.w.communicate< Edge, Vertex >( level );
   }

   this->timingTree_->start( "Macro-Vertex" );

   for ( const auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.u.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            vertexdof::macrovertex::projectNormal3D< real_t >(level, vertex, storage_, normal_function_, dst.u.getVertexDataID(), dst.v.getVertexDataID(), dst.w.getVertexDataID() );
         }
         else
         {
            vertexdof::macrovertex::projectNormal2D< real_t >(level, vertex, storage_, normal_function_, dst.u.getVertexDataID(), dst.v.getVertexDataID() );
         }
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );

   this->timingTree_->start( "Macro-Edge" );

   if ( level >= 1 )
   {
     for ( const auto & it : storage_->getEdges())
     {
       Edge & edge = *it.second;

       const DoFType edgeBC = dst.u.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag());
       if ( testFlag( edgeBC, flag ))
       {
          if ( storage_->hasGlobalCells() )
          {
             vertexdof::macroedge::projectNormal3D< real_t >(level, edge, storage_, normal_function_, dst.u.getEdgeDataID(), dst.v.getEdgeDataID(), dst.w.getEdgeDataID() );
          }
          else
          {
             vertexdof::macroedge::projectNormal2D< real_t >(level, edge, storage_, normal_function_, dst.u.getEdgeDataID(), dst.v.getEdgeDataID() );
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

         const DoFType faceBC = dst.u.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroface::projectNormal3D< real_t >(level, face, storage_, normal_function_, dst.u.getFaceDataID(), dst.v.getFaceDataID(), dst.w.getFaceDataID() );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   this->stopTiming( "Apply" );
}

} // namespace hyteg
