#pragma once

#include <array>

namespace hhg {
namespace P1 {

template<uint_t N>
inline real_t getSupremumNorm(std::array<P1Function<real_t>*, N> functions, uint_t level) {
  WALBERLA_ASSERT(N >= 1, "getSupremumNorm requires at least one function as argument");

  real_t maxValue = functions[0]->getMaxValue(level);

  for (uint_t k = 1; k < functions.size(); ++k) {
    maxValue = std::max(maxValue, functions[k]->getMaxValue(level));
  }

  return maxValue;
}

template<uint_t N>
inline real_t getApproximateEuclideanNorm(std::array<P1Function<real_t>*, N> functions, uint_t level) {
  return std::sqrt(real_c(N)) * getSupremumNorm<N>(functions, level);
}

template<uint_t N>
inline real_t projectNormal(const std::shared_ptr< PrimitiveStorage > & storage, std::array<P1Function<real_t>*, N> velocity, std::array<P1Function<real_t>*, N> normals, uint_t level, DoFType flag) {

  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Vertex>> velocityVertexIDs;
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Edge>>   velocityEdgeIDs;

  for (auto& v : velocity) {
    velocityVertexIDs.push_back(v->getVertexDataID());
    velocityEdgeIDs.push_back(v->getEdgeDataID());
  }

  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Vertex>> normalVertexIDs;
  std::vector<PrimitiveDataID<FunctionMemory< real_t >, Edge>>   normalEdgeIDs;

  for (auto& n : normals) {
    normalVertexIDs.push_back(n->getVertexDataID());
    normalEdgeIDs.push_back(n->getEdgeDataID());
  }

  for (auto& it : storage->getVertices()) {
    Vertex& vertex = *it.second;

    if (testFlag(vertex.getDoFType(), flag)) {
      P1Vertex::projectNormal(vertex, velocityVertexIDs, normalVertexIDs, level);
    }
  }

  for (auto& v : velocity) {
    v->getCommunicator(level)->template startCommunication<Vertex, Edge>();
  }

  for (auto& it : storage->getEdges()) {
    Edge& edge = *it.second;

    if (testFlag(edge.getDoFType(), flag)) {
      P1Edge::projectNormal<real_t>(level, edge, velocityEdgeIDs, normalEdgeIDs);
    }
  }

  for (auto& v : velocity) {
    v->getCommunicator(level)->template endCommunication<Vertex, Edge>();
  }

  for (auto& v : velocity) {
    v->getCommunicator(level)->template startCommunication<Edge, Face>();
    v->getCommunicator(level)->template endCommunication<Edge, Face>();
  }
}

}
}