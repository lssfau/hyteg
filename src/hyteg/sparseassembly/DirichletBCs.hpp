/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P2StokesFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

// ==================
//  Scalar Functions
// ==================
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

inline void applyDirichletBC( const EdgeDoFFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( auto& it : numerator.getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = numerator.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, DirichletBoundary ) )
      {
         edgedof::macroedge::applyDirichletBC( level, edge, mat, numerator.getEdgeDataID() );
      }
   }

   for ( auto& it : numerator.getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = numerator.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, DirichletBoundary ) )
      {
         edgedof::macroface::applyDirichletBC( level, face, mat, numerator.getFaceDataID() );
      }
   }
}

inline void applyDirichletBC( const P2Function< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   applyDirichletBC( numerator.getVertexDoFFunction(), mat, level );
   applyDirichletBC( numerator.getEdgeDoFFunction(), mat, level );
}

// ==================
//  Vector Functions
// ==================
inline void applyDirichletBC( const P1VectorFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( uint_t k = 0; k < numerator.getDimension(); k++ )
   {
      applyDirichletBC( numerator[k], mat, level );
   }
}

inline void applyDirichletBC( const P2VectorFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( uint_t k = 0; k < numerator.getDimension(); k++ )
   {
      applyDirichletBC( numerator[k], mat, level );
   }
}

// ==================
//  Stokes Functions
// ==================
inline void applyDirichletBC( const P1StokesFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( uint_t k = 0; k < numerator.uvw().getDimension(); ++k )
   {
      applyDirichletBC( numerator.uvw()[k], mat, level );
   }
}

inline void applyDirichletBC( const P2P1TaylorHoodFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( uint_t k = 0; k < numerator.uvw().getDimension(); ++k )
   {
      applyDirichletBC( numerator.uvw()[k], mat, level );
   }
}

inline void applyDirichletBC( const P2P2StokesFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( uint_t k = 0; k < numerator.uvw().getDimension(); ++k )
   {
      applyDirichletBC( numerator.uvw()[k], mat, level );
   }
}

} // namespace hyteg
