#include "P1Function.hpp"
#include "P1DataHandling.hpp"
#include "p1vertex.hpp"
#include "p1edge.hpp"
#include "p1face.hpp"

namespace hhg {

P1Function::P1Function(const std::string& name, PrimitiveStorage& storage, uint_t minLevel, uint_t maxLevel)
    : Function(name, storage, minLevel, maxLevel)
{
    FaceP1FunctionMemoryDataHandling faceP1FunctionMemoryDataHandling(minLevel, maxLevel);
    EdgeP1FunctionMemoryDataHandling edgeP1FunctionMemoryDataHandling(minLevel, maxLevel);
    VertexP1FunctionMemoryDataHandling vertexP1FunctionMemoryDataHandling(minLevel, maxLevel);
    faceDataID_ = storage.addFaceData(faceP1FunctionMemoryDataHandling, name);
    edgeDataID_ = storage.addEdgeData(edgeP1FunctionMemoryDataHandling, name);
    vertexDataID_ = storage.addVertexData(vertexP1FunctionMemoryDataHandling, name);
}

P1Function::~P1Function()
{
    //TODO implement!
}

void P1Function::interpolate(std::function<real_t(const hhg::Point3D&)>& expr, uint_t level, DoFType flag = All)
{
    for (auto& it : storage_.getVertices()) {
        Vertex& vertex = *it.second;

        if (testFlag(vertex.type, flag)) {
            P1Vertex::interpolate(vertex, vertexDataID_, expr, level);
        }
    }

    comm[level].startCommunication<Vertex, Edge>();

    for (auto& it : storage_.getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.type, flag)) {
            P1Edge::interpolate(edge, edgeDataID_, expr, level);
        }
    }

    comm[level].endCommunication<Vertex, Edge>();
    comm[level].startCommunication<Edge, Face>();

    for (auto& it : storage_.getFaces()) {
        Face& face = *it.second;

        if (testFlag(face.type, flag)) {
            P1Face::interpolate(face, faceDataID_, expr, level);
        }
    }

    comm[level].endCommunication<Edge, Face>();
}

void P1Function::assign(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag = All)
{
    // Collect all source IDs in a vector
    std::vector<const PrimitiveDataID<VertexP1FunctionMemory, Vertex>> srcVertexIDs(functions.size());
    std::vector<const PrimitiveDataID<EdgeP1FunctionMemory, Edge>> srcEdgeIDs(functions.size());
    std::vector<const PrimitiveDataID<FaceP1FunctionMemory, Face>> srcFaceIDs(functions.size());

    for (auto& function : functions)
    {
        srcVertexIDs.push_back(function->vertexDataID_);
        srcEdgeIDs.push_back(function->edgeDataID_);
        srcFaceIDs.push_back(function->faceDataID_);
    }

    for (auto& it : storage_.getVertices()) {
        Vertex& vertex = *it.second;

        if (testFlag(vertex.type, flag)) {
            P1Vertex::assign(vertex, scalars, srcVertexIDs, vertexDataID_, level);
        }
    }

    comm[level].startCommunication<Vertex, Edge>();

    for (auto& it : storage_.getEdges()) {
        Edge& edge = *it.second;

        if (testFlag(edge.type, flag)) {
            P1Edge::assign(edge, scalars, srcEdgeIDs, edgeDataID_, level);
        }
    }

    comm[level].endCommunication<Vertex, Edge>();
    comm[level].startCommunication<Edge, Face>();

    for (auto& it : storage_.getFaces()) {
        Face& face = *it.second;

        if (testFlag(face.type, flag)) {
            P1Face::assign(level, face, scalars, srcFaceIDs, faceDataID_);
        }
    }

    comm[level].endCommunication<Edge, Face>();
}

void P1Function::add(const std::vector<walberla::real_t> scalars, const std::vector<P1Function*> functions, size_t level, DoFType flag = All)
{

}

real_t P1Function::dot(P1Function& rhs, size_t level, DoFType flag = All)
{

}

void P1Function::prolongate(size_t level, DoFType flag = All){

}

void P1Function::prolongateQuadratic(size_t level, DoFType flag = All){

}

void P1Function::restrict(size_t level, DoFType flag = All)
{

}
}
