#pragma once

#include "tinyhhg_core/p1functionspace/P1Vertex.hpp"
#include "tinyhhg_core/p1functionspace/P1Edge.hpp"
#include "tinyhhg_core/p1functionspace/P1Face.hpp"

namespace hhg {
namespace petsc {

inline void createVectorFromFunction(P1Function<PetscScalar> &function,
                                     P1Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  for (auto &it : function.getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      vertexdof::macrovertex::createVectorFromFunction<PetscScalar>(vertex, function.getVertexDataID(), numerator.getVertexDataID(), vec, level);
    }
  }

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      P1Edge::createVectorFromFunction<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      vertexdof::macroface::createVectorFromFunction<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }
}

inline void createFunctionFromVector(P1Function<PetscScalar> &function,
                                     P1Function<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  for (auto &it : function.getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      vertexdof::macrovertex::createFunctionFromVector<PetscScalar>(vertex, function.getVertexDataID(), numerator.getVertexDataID(), vec, level);
    }
  }

  function.getCommunicator(level)->template startCommunication<Vertex, Edge>();
  function.getCommunicator(level)->template endCommunication<Vertex, Edge>();

  for (auto &it : function.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      P1Edge::createFunctionFromVector<PetscScalar>(level, edge, function.getEdgeDataID(), numerator.getEdgeDataID(), vec);
    }
  }

  function.getCommunicator(level)->template startCommunication<Edge, Face>();
  function.getCommunicator(level)->template endCommunication<Edge, Face>();

  for (auto &it : function.getStorage()->getFaces()) {
    Face &face = *it.second;

    if (testFlag(face.type, flag)) {
      vertexdof::macroface::createFunctionFromVector<PetscScalar>(level, face, function.getFaceDataID(), numerator.getFaceDataID(), vec);
    }
  }
}

inline void applyDirichletBC(P1Function<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  for (auto &it : numerator.getStorage()->getVertices()) {
    Vertex &vertex = *it.second;

    if (testFlag(vertex.getDoFType(), DirichletBoundary)) {
      vertexdof::macrovertex::applyDirichletBC(vertex, mat, level, numerator.getVertexDataID());
    }
  }

  for (auto &it : numerator.getStorage()->getEdges()) {
    Edge &edge = *it.second;

    if (testFlag(edge.getDoFType(), DirichletBoundary)) {
      P1Edge::applyDirichletBC(level, edge, mat, numerator.getEdgeDataID());
    }
  }

}

template<class OperatorType>
inline void createMatrix(OperatorType& opr, P1Function< PetscInt > & src, P1Function< PetscInt > & dst, Mat& mat, size_t level, DoFType flag)
{
  for (auto& it : opr.getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      vertexdof::macrovertex::saveOperator(vertex, opr.getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat, level);
    }
  }

  for (auto& it : opr.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      P1Edge::saveOperator(level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
    }
  }

  for (auto& it : opr.getStorage()->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      vertexdof::macroface::saveOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
    }
  }
}

}
}
