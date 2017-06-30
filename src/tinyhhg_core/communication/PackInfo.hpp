
#pragma once

namespace hhg {
/// namespace containing function for communication between primitives
namespace communication {
// WIP
///
/// /brief Abstract class for pack and unpack functions of primitives
///
/// This abstract class provides the necessary functions to handle communication between the different primitives
/// (Vertex, Edge, Face). For each pair of primitives which need to communicate,
/// there must be a pack, unpack and a communicateLocal function
///
class PackInfo {
public:

    virtual ~PackInfo() {}

    //virtual void packDataEdgeToFace        ( const Edge * sender,   const PrimitiveID & receiver, walberla::mpi::SendBuffer & buffer ) = 0;
    //virtual void unpackDataEdgeToFace      (       Face * receiver, const PrimitiveID & sender  , walberla::mpi::RecvBuffer & buffer ) = 0;

    //virtual void communicateLocalEdgeToFace( const Edge * sender, Face * receiver ) = 0;

};

} // namespace communication
} // namespace hhg
