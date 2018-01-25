#pragma once

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/StencilMemory.hpp"

#include "tinyhhg_core/levelinfo.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/polynomial/Polynomial2D.hpp"

#include <string>


namespace hhg
{

/////////////////////
// Function memory //
/////////////////////

inline uint_t P1VertexFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  return levelinfo::num_microvertices_per_vertex( level ) + numDependencies;
}

inline uint_t P1EdgeFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  size_t num_dofs_per_edge = levelinfo::num_microvertices_per_edge( level );
  return num_dofs_per_edge + numDependencies * ( num_dofs_per_edge - 1 );
}

inline uint_t P1FaceFunctionMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  return levelinfo::num_microvertices_per_face(level);
}

template< typename ValueType >
using VertexP1FunctionMemory = FunctionMemory< ValueType >;

template< typename ValueType >
using EdgeP1FunctionMemory = FunctionMemory< ValueType >;

template< typename ValueType >
using FaceP1FunctionMemory = FunctionMemory< ValueType >;


////////////////////
// Stencil memory //
////////////////////

template< typename ValueType >
using VertexP1StencilMemory = StencilMemory< ValueType >;

template< typename ValueType >
using EdgeP1StencilMemory = StencilMemory< ValueType >;

template< typename ValueType >
using FaceP1StencilMemory = StencilMemory< ValueType >;

inline uint_t P1VertexStencilMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  return levelinfo::num_microvertices_per_vertex( level ) + numDependencies;
}

inline uint_t P1EdgeStencilMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( numDependencies );
  return 7;
}

inline uint_t P1FaceStencilMemorySize( const uint_t & level, const uint_t & numDependencies )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( numDependencies );
  return 7;
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
    std::shared_ptr<GeneralPolynomial2D> horiPoly;
    std::shared_ptr<GeneralPolynomial2D> vertPoly;
    std::shared_ptr<GeneralPolynomial2D> diagPoly;
  };

  std::map<uint_t, FacePolynomials> polynomials_;

  FaceP1PolynomialMemory(uint_t maxDegree) {
    addDegree(maxDegree);
  }

  inline FacePolynomials& addDegree(uint_t degree)
  {
    if (polynomials_.count(degree)>0) {
      WALBERLA_LOG_WARNING("Degree already exists.")
    }

    FacePolynomials& tmp = polynomials_[degree];
    tmp.horiPoly = std::make_shared<GeneralPolynomial2D>(degree);
    tmp.vertPoly = std::make_shared<GeneralPolynomial2D>(degree);
    tmp.diagPoly = std::make_shared<GeneralPolynomial2D>(degree);
    return tmp;
  }

  GeneralPolynomial2D& getHoriPolynomial(uint_t maxDegree) {
    return *polynomials_[maxDegree].horiPoly;
  }

  GeneralPolynomial2D& getVertPolynomial(uint_t maxDegree) {
    return *polynomials_[maxDegree].vertPoly;
  }

  GeneralPolynomial2D& getDiagPolynomial(uint_t maxDegree) {
    return *polynomials_[maxDegree].diagPoly;
  }

};


}
