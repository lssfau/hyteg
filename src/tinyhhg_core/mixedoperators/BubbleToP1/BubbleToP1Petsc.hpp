#pragma once

namespace hhg {
namespace petsc {

template<class OperatorType>
inline void createMatrix(OperatorType& opr, BubbleFunction< PetscInt > & src, P1Function< PetscInt > & dst, Mat& mat, uint_t level, DoFType flag)
{
  src.getCommunicator(level)->startCommunication<Face, Edge>();
  src.getCommunicator(level)->endCommunication<Face, Edge>();

  src.getCommunicator(level)->startCommunication<Edge, Vertex>();
  src.getCommunicator(level)->endCommunication<Edge, Vertex>();

  for (auto& it : src.getStorage()->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag))
    {
      BubbleToP1Vertex::saveOperator(vertex, opr.getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat, level);
    }
  }

  for (auto& it : src.getStorage()->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag))
    {
      BubbleToP1Edge::saveOperator(level, edge, opr.getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat);
    }
  }

  for (auto& it : src.getStorage()->getFaces()) {
    Face& face = *it.second;

    if (testFlag(face.type, flag))
    {
      BubbleToP1Face::saveOperator(level, face, opr.getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat);
    }
  }
}

}
}