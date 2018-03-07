
#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroFace.hpp"

namespace hhg {
namespace EdgeDoF {

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

    if (testFlag(edge.getDoFType(), flag)) {
      edgedof::macroedge::createVectorFromFunction<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      edgedof::macroface::createVectorFromFunction<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }
}

inline void createFunctionFromVector(EdgeDoFFunction<PetscScalar> &function,
                                     EdgeDoFFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  function.getCommunicator(level)->template startCommunication<Vertex, Edge>();
  function.getCommunicator(level)->template endCommunication<Vertex, Edge>();

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      edgedof::macroedge::createFunctionFromVector<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  function.getCommunicator(level)->template startCommunication<Edge, Face>();
  function.getCommunicator(level)->template endCommunication<Edge, Face>();

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      edgedof::macroface::createFunctionFromVector<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }
}

inline void applyDirichletBC(EdgeDoFFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {

  for (auto &it : numerator.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), DirichletBoundary)) {
      edgedof::macroedge::applyDirichletBC(level, edge, mat, numerator.getEdgeDataID());
    }
  }

}


template<uint_t Level>
inline void saveEdgeOperatorTmpl( const Edge & edge,
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

SPECIALIZE(void, saveEdgeOperatorTmpl, saveEdgeOperator);

template<uint_t Level>
inline void saveFaceOperatorTmpl( const Face & face,
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

SPECIALIZE(void, saveFaceOperatorTmpl, saveFaceOperator);

template<class OperatorType>
inline void createMatrix(OperatorType& opr, EdgeDoFFunction< PetscInt > & src, EdgeDoFFunction< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  for (auto& it : opr.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      saveEdgeOperator(level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
    }
  }

  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      saveFaceOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
    }
  }
}

#endif

}

}