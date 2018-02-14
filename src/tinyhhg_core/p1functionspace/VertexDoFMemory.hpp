#pragma once

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/types/matrix.hpp"

#include <string>


namespace hhg {

/////////////////////
// Function memory //
/////////////////////

inline uint_t vertexDoFMacroVertexFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  return levelinfo::num_microvertices_per_vertex( level ) + primitive.getNumNeighborEdges();
}

inline uint_t vertexDoFMacroEdgeFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  size_t num_dofs_per_edge = levelinfo::num_microvertices_per_edge( level );
  return num_dofs_per_edge + primitive.getNumNeighborFaces() * ( num_dofs_per_edge - 1 );
}

inline uint_t vertexDoFMacroFaceFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  return levelinfo::num_microvertices_per_face( level );
}

inline uint_t vertexDoFMacroCellFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  return levelinfo::num_microvertices_per_cell( level );
}


////////////////////
// Stencil memory //
////////////////////

inline uint_t vertexDoFMacroVertexStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  return levelinfo::num_microvertices_per_vertex( level ) + primitive.getNumNeighborEdges();
}

inline uint_t vertexDoFMacroEdgeStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 7;
}

inline uint_t vertexDoFMacroFaceStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 7;
}

inline uint_t vertexDoFMacroCellStencilMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 15;
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

} // namespace hhg
