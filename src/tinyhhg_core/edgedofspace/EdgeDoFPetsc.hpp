
#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroFace.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

namespace hhg {
namespace edgedof {

using walberla::real_t;
using walberla::uint_t;

#ifdef HHG_BUILD_WITH_PETSC

inline void createVectorFromFunction(EdgeDoFFunction<PetscScalar> &function,
                                     EdgeDoFFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      macroedge::createVectorFromFunction<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
      macroface::createVectorFromFunction<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag)) {
      macrocell::createVectorFromFunction<PetscScalar>(level, cell, function.getCellDataID(), numerator.getCellDataID(), vec);
    }
  }
}

inline void createFunctionFromVector(EdgeDoFFunction<PetscScalar> &function,
                                     EdgeDoFFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  function.startCommunication<Vertex, Edge>( level );
  function.endCommunication<Vertex, Edge>( level );

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag)) {
      edgedof::macroedge::createFunctionFromVector<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  function.startCommunication<Edge, Face>( level );
  function.endCommunication<Edge, Face>( level );

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag)) {
      edgedof::macroface::createFunctionFromVector<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getCells()) {
    Cell & cell = *it.second;

    const DoFType cellBC = function.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag)) {
      edgedof::macrocell::createFunctionFromVector<PetscScalar>(level, cell, function.getCellDataID(), numerator.getCellDataID(), vec);
    }
  }
}

inline void applyDirichletBC(EdgeDoFFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {

  for (auto &it : numerator.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    const DoFType edgeBC = numerator.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, DirichletBoundary)) {
      edgedof::macroedge::applyDirichletBC(level, edge, mat, numerator.getEdgeDataID());
    }
  }

}



inline void saveEdgeOperator( const uint_t & Level, const Edge & edge,
                              const PrimitiveDataID< StencilMemory< real_t >, Edge>    & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Edge> & dstId,
                              Mat & mat )
{
  using namespace hhg::edgedof::macroedge;
  size_t rowsize = levelinfo::num_microedges_per_edge(Level);

  real_t * opr_data = edge.getData(operatorId)->getPointer( Level );
  PetscInt * src      = edge.getData(srcId)->getPointer( Level );
  PetscInt * dst      = edge.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetscInt dstInt;

  for(uint_t i = 0; i < rowsize; ++i){
    dstInt = dst[indexFromHorizontalEdge( Level, i, stencilDirection::EDGE_HO_C )];

    for(uint_t k = 0; k < neighborsOnEdgeFromHorizontalEdge.size(); ++k){
      srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnEdgeFromHorizontalEdge[k] )];
      MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnEdgeFromHorizontalEdge[k])], INSERT_VALUES);
    }
    for(uint_t k = 0; k < neighborsOnSouthFaceFromHorizontalEdge.size(); ++k){
      srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnSouthFaceFromHorizontalEdge[k] )];
      MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnSouthFaceFromHorizontalEdge[k])], INSERT_VALUES);
    }
    if(edge.getNumNeighborFaces() == 2){
      for(uint_t k = 0; k < neighborsOnNorthFaceFromHorizontalEdge.size(); ++k){
        srcInt = src[indexFromHorizontalEdge( Level, i, neighborsOnNorthFaceFromHorizontalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[hhg::edgedof::stencilIndexFromHorizontalEdge(neighborsOnNorthFaceFromHorizontalEdge[k])], INSERT_VALUES);
      }
    }
  }
}



inline void saveFaceOperator( const uint_t & Level, const Face & face,
                              const PrimitiveDataID< StencilMemory< real_t >, Face>    & operatorId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & srcId,
                              const PrimitiveDataID< FunctionMemory< PetscInt >, Face> & dstId,
                              Mat & mat )
{
  real_t * opr_data = face.getData(operatorId)->getPointer( Level );
  PetscInt * src      = face.getData(srcId)->getPointer( Level );
  PetscInt * dst      = face.getData(dstId)->getPointer( Level );

  PetscInt srcInt;
  PetscInt dstInt;

  using namespace edgedof::macroface;

  for ( const auto & it : hhg::edgedof::macroface::Iterator( Level, 0 ) )
  {
    if( it.row() != 0) {
      dstInt = dst[indexFromHorizontalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_HO_C )];
      for(uint_t k = 0; k < neighborsFromHorizontalEdge.size(); ++k){
        srcInt = src[indexFromHorizontalEdge( Level, it.col(), it.row(), neighborsFromHorizontalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[edgedof::stencilIndexFromHorizontalEdge(neighborsFromHorizontalEdge[k])], INSERT_VALUES);
      }
    }
    if( it.col() + it.row() != (hhg::levelinfo::num_microedges_per_edge( Level ) - 1)) {
      dstInt = dst[indexFromDiagonalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_DI_C )];
      for(uint_t k = 0; k < neighborsFromDiagonalEdge.size(); ++k){
        srcInt = src[indexFromDiagonalEdge( Level, it.col(), it.row(), neighborsFromDiagonalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[edgedof::stencilIndexFromDiagonalEdge(neighborsFromDiagonalEdge[k])], INSERT_VALUES);
      }
    }
    if( it.col() != 0) {
      dstInt = dst[indexFromVerticalEdge( Level, it.col(), it.row(), stencilDirection::EDGE_VE_C )];
      for(uint_t k = 0; k < neighborsFromVerticalEdge.size(); ++k){
        srcInt = src[indexFromVerticalEdge( Level, it.col(), it.row(), neighborsFromVerticalEdge[k] )];
        MatSetValues(mat, 1, &dstInt, 1, &srcInt, &opr_data[edgedof::stencilIndexFromVerticalEdge(neighborsFromVerticalEdge[k])], INSERT_VALUES);
      }
    }
  }
}


template<class OperatorType>
inline void createMatrix(OperatorType& opr, EdgeDoFFunction< PetscInt > & src, EdgeDoFFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  for (auto& it : opr.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag))
    {
      saveEdgeOperator(level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
    }
  }

  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      saveFaceOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
    }
  }
}

#endif

}

}