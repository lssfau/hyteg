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

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {
namespace petsc {

// ============
//  P1Function
// ============
inline void createVectorFromFunction( const P1Function< PetscReal >&        function,
                                      const P1Function< idx_t >&            numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   for ( auto& it : function.getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::createVectorFromFunction< PetscReal >(
             vertex, function.getVertexDataID(), numerator.getVertexDataID(), vec, level );
      }
   }

   for ( auto& it : function.getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::createVectorFromFunction< PetscReal >(
             level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec );
      }
   }

   for ( auto& it : function.getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         vertexdof::macroface::createVectorFromFunction< PetscReal >(
             level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec );
      }
   }

   for ( auto& it : function.getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         vertexdof::macrocell::createVectorFromFunction< PetscReal >(
             level, cell, function.getCellDataID(), numerator.getCellDataID(), vec );
      }
   }
}

inline void createFunctionFromVector( const P1Function< PetscReal >&        function,
                                      const P1Function< idx_t >&            numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   for ( auto& it : function.getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::createFunctionFromVector< PetscReal >(
             vertex, function.getVertexDataID(), numerator.getVertexDataID(), vec, level );
      }
   }

   function.startCommunication< Vertex, Edge >( level );
   function.endCommunication< Vertex, Edge >( level );

   for ( auto& it : function.getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::createFunctionFromVector< PetscReal >(
             level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec );
      }
   }

   function.startCommunication< Edge, Face >( level );
   function.endCommunication< Edge, Face >( level );

   for ( auto& it : function.getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         vertexdof::macroface::createFunctionFromVector< PetscReal >(
             level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec );
      }
   }

   function.startCommunication< Face, Cell >( level );
   function.endCommunication< Face, Cell >( level );

   for ( auto& it : function.getStorage()->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         vertexdof::macrocell::createFunctionFromVector< PetscReal >(
             level, cell, function.getCellDataID(), numerator.getCellDataID(), vec );
      }
   }
}

inline void applyDirichletBC( const P1Function< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( auto& it : numerator.getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = numerator.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, DirichletBoundary ) )
      {
         vertexdof::macrovertex::applyDirichletBC( vertex, mat, level, numerator.getVertexDataID() );
      }
   }

   for ( auto& it : numerator.getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = numerator.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, DirichletBoundary ) )
      {
         vertexdof::macroedge::applyDirichletBC( level, edge, mat, numerator.getEdgeDataID() );
      }
   }

   for ( auto& it : numerator.getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = numerator.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, DirichletBoundary ) )
      {
         vertexdof::macroface::applyDirichletBC( level, face, mat, numerator.getFaceDataID() );
      }
   }
}

// ==================
//  P1VectorFunction
// ==================
inline void createVectorFromFunction( const P1VectorFunction< PetscReal >&  function,
                                      const P1VectorFunction< idx_t >&      numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   for ( uint_t k = 0; k < function.getDimension(); k++ )
   {
      createVectorFromFunction( function[k], numerator[k], vec, level, flag );
   }
}

inline void createFunctionFromVector( const P1VectorFunction< PetscReal >&  function,
                                      const P1VectorFunction< idx_t >&      numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   for ( uint_t k = 0; k < function.getDimension(); k++ )
   {
      createFunctionFromVector( function[k], numerator[k], vec, level, flag );
   }
}

inline void applyDirichletBC( const P1VectorFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( uint_t k = 0; k < numerator.getDimension(); k++ )
   {
      applyDirichletBC( numerator[k], mat, level );
   }
}

// =============
//  P1Operators
// =============
inline void saveIdentityOperator( const P1Function< idx_t >&                  dst,
                                  const std::shared_ptr< SparseMatrixProxy >& mat,
                                  size_t                                      level,
                                  DoFType                                     flag )
{
   const auto storage = dst.getStorage();

   for ( auto& it : storage->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::saveIdentityOperator( vertex, dst.getVertexDataID(), mat, level );
      }
   }

   for ( auto& it : storage->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::saveIdentityOperator( level, edge, dst.getEdgeDataID(), mat );
      }
   }

   if ( level >= 2 )
   {
      for ( auto& it : storage->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::saveIdentityOperator( level, face, dst.getFaceDataID(), mat );
         }
      }

      for ( auto& it : storage->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            vertexdof::macrocell::saveIdentityOperator( level, cell, dst.getCellDataID(), mat );
         }
      }
   }
}

} // namespace petsc
} // namespace hyteg

#endif
