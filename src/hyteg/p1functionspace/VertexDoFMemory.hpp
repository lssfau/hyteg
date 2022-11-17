/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "core/Abort.h"

#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"

#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/polynomial/Polynomial.hpp"

#include <string>


namespace hyteg {

/////////////////////
// Function memory //
/////////////////////

inline uint_t vertexDoFMacroVertexFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  return levelinfo::num_microvertices_per_vertex( level ) + primitive.getNumNeighborEdges();
}

inline uint_t vertexDoFMacroEdgeFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  const uint_t numDoFsPerEdge = levelinfo::num_microvertices_per_edge( level );
  const uint_t numDoFsPerFaceGhostLayer = levelinfo::num_microvertices_per_edge( level ) - 1;
  const uint_t numDoFsPerCellGhostLayer = levelinfo::num_microvertices_per_edge( level ) - 2;
  return numDoFsPerEdge
         + primitive.getNumNeighborFaces() * numDoFsPerFaceGhostLayer
         + primitive.getNumNeighborCells() * numDoFsPerCellGhostLayer;
}

inline uint_t vertexDoFMacroFaceFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  const uint_t width = (1 << level) + 1;
  return levelinfo::num_microvertices_per_face( level ) + primitive.getNumNeighborCells() * levelinfo::num_microvertices_per_face_from_width( width - 1 );
}

inline uint_t vertexDoFMacroCellFunctionMemorySize( const uint_t & level, const Cell & primitive )
{
  WALBERLA_UNUSED( primitive );
  return levelinfo::num_microvertices_per_cell( level );
}

inline unsigned long long vertexDoFLocalFunctionMemorySize( const uint_t & level, const std::shared_ptr< PrimitiveStorage > & storage )
{
   unsigned long long mem = 0;

   for ( const auto & it : storage->getVertices() )
   {
      mem += vertexDoFMacroVertexFunctionMemorySize( level, *it.second );
   }

   for ( const auto & it : storage->getEdges() )
   {
      mem += vertexDoFMacroEdgeFunctionMemorySize( level, *it.second );
   }

   for ( const auto & it : storage->getFaces() )
   {
      mem += vertexDoFMacroFaceFunctionMemorySize( level, *it.second );
   }

   for ( const auto & it : storage->getCells() )
   {
      mem += vertexDoFMacroCellFunctionMemorySize( level, *it.second );
   }

   return mem;
}

inline unsigned long long vertexDoFGlobalFunctionMemorySize( const uint_t & level, const std::shared_ptr< PrimitiveStorage > & storage )
{
   const auto memLocal = vertexDoFLocalFunctionMemorySize( level, storage );
   const auto memGlobal = walberla::mpi::allReduce( memLocal, walberla::mpi::SUM );
   return memGlobal;
}

////////////////////
// Stencil memory //
////////////////////

// constant stencils

inline uint_t vertexDoFMacroVertexStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  return levelinfo::num_microvertices_per_vertex( level ) + primitive.getNumNeighborEdges();
}

inline uint_t vertexDoFMacroEdgeStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return 3 + 2 * primitive.getNumNeighborFaces() + primitive.getNumNeighborCells();
}

inline uint_t vertexDoFMacroFaceStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 27;
}

inline uint_t vertexDoFMacroCellStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 27;
}

// pointwise stencils

inline uint_t vertexDoFMacroVertexPointwiseStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  return levelinfo::num_microvertices_per_vertex( level ) + primitive.getNumNeighborEdges();
}

inline uint_t vertexDoFMacroEdgePointwiseStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  if ( primitive.getNumNeighborCells() == 0 )
  {
    const uint_t numDoFs = levelinfo::num_microvertices_per_edge( level ) - 2;
    return numDoFs * (3 + 2 * primitive.getNumNeighborFaces());
  }
  WALBERLA_ABORT( "Point-wise stencils not implemented for 3D." );
}

inline uint_t vertexDoFMacroFacePointwiseStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  if ( primitive.getNumNeighborCells() == 0 )
  {
    const uint_t numDoFs = levelinfo::num_microvertices_per_face( level ) -
      3 * (levelinfo::num_microvertices_per_edge( level ) - 1);
    return numDoFs * 7;
  }
  WALBERLA_ABORT( "Point-wise stencils not implemented for 3D." );
}

class VertexP1LocalMatrixMemory
{
public:

  std::map<size_t, std::vector<Matrix3r>> data;
  size_t num_deps_;

  inline std::vector<Matrix3r>& addlevel(size_t level, size_t num_deps)
  {
    if (data.count(level)>0)
    WALBERLA_LOG_WARNING("Level already exists.")
    else
    {
      this->num_deps_ = num_deps;
      data[level].resize(num_deps);
    }
    return data[level];
  }

  Matrix3r& getGrayMatrix(uint_t level, uint_t idx) {
    return data[level][idx];
  }

  inline size_t getSize(size_t level)
  {
    return num_deps_;
  }

};

class EdgeP1LocalMatrixMemory
{
public:

  std::map<size_t, std::array<Matrix3r, 2>> dataGray;
  std::map<size_t, std::array<Matrix3r, 2>> dataBlue;

  inline std::array<Matrix3r, 2>& addlevel(size_t level)
  {
    if (dataGray.count(level)>0) {
      WALBERLA_LOG_WARNING("Level already exists.")
    }
    return dataGray[level];
  }

  Matrix3r& getGrayMatrix(uint_t level, uint_t idx) {
    return dataGray[level][idx];
  }

  Matrix3r& getBlueMatrix(uint_t level, uint_t idx) {
    return dataBlue[level][idx];
  }

};

class FaceP1LocalMatrixMemory
{
public:

  std::map<size_t, std::array<Matrix3r, 2>> data;

  inline std::array<Matrix3r, 2>& addlevel(size_t level)
  {
    if (data.count(level)>0) {
      WALBERLA_LOG_WARNING("Level already exists.")
    }
    return data[level];
  }

  Matrix3r& getGrayMatrix(uint_t level) {
    return data[level][0];
  }

  Matrix3r& getBlueMatrix(uint_t level) {
    return data[level][1];
  }

};

class FaceP1PolynomialMemory
{
public:

  struct FacePolynomials {
    std::shared_ptr<GeneralPolynomial2D> polyS;
    std::shared_ptr<GeneralPolynomial2D> polySE;
    std::shared_ptr<GeneralPolynomial2D> polyW;
    std::shared_ptr<GeneralPolynomial2D> polyC;
    std::shared_ptr<GeneralPolynomial2D> polyE;
    std::shared_ptr<GeneralPolynomial2D> polyNW;
    std::shared_ptr<GeneralPolynomial2D> polyN;
  };

  std::map<uint_t, FacePolynomials> polynomials_;

  FaceP1PolynomialMemory() {}

  inline FacePolynomials& addDegree(uint_t degree)
  {
    if (polynomialDegreeExists(degree)) {
      WALBERLA_LOG_WARNING("Degree already exists.");
    }

    FacePolynomials& tmp = polynomials_[degree];
    tmp.polyS = std::make_shared<GeneralPolynomial2D>(degree);
    tmp.polySE = std::make_shared<GeneralPolynomial2D>(degree);
    tmp.polyW = std::make_shared<GeneralPolynomial2D>(degree);
    tmp.polyC = std::make_shared<GeneralPolynomial2D>(degree);

    // TODO: only required in asymmetric case
    tmp.polyE = std::make_shared<GeneralPolynomial2D>(degree);
    tmp.polyNW = std::make_shared<GeneralPolynomial2D>(degree);
    tmp.polyN = std::make_shared<GeneralPolynomial2D>(degree);
    return tmp;
  }

  bool polynomialDegreeExists(uint_t degree) {
    return polynomials_.count(degree)>0;
  }

  GeneralPolynomial2D& getPolynomialS(uint_t maxDegree) {
    return *polynomials_[maxDegree].polyS;
  }

  GeneralPolynomial2D& getPolynomialSE(uint_t maxDegree) {
    return *polynomials_[maxDegree].polySE;
  }

  GeneralPolynomial2D& getPolynomialW(uint_t maxDegree) {
    return *polynomials_[maxDegree].polyW;
  }

  GeneralPolynomial2D& getPolynomialC(uint_t maxDegree) {
    return *polynomials_[maxDegree].polyC;
  }

  GeneralPolynomial2D& getPolynomialE(uint_t maxDegree) {
    return *polynomials_[maxDegree].polyE;
  }

  GeneralPolynomial2D& getPolynomialNW(uint_t maxDegree) {
    return *polynomials_[maxDegree].polyNW;
  }

  GeneralPolynomial2D& getPolynomialN(uint_t maxDegree) {
    return *polynomials_[maxDegree].polyN;
  }

};


} // namespace hyteg
