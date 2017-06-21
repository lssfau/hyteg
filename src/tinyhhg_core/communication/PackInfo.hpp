
#pragma once

namespace hhg {

// WIP

class PackInfo
{
public:

  virtual ~PackInfo() {}

  virtual void turnAfterPacking( bool lol );

  virtual void packDataEdgeToFace        ( const Edge * sender,   const PrimitiveID & receiver, walberla::mpi::SendBuffer & buffer ) = 0;
  virtual void unpackDataEdgeToFace      (       Face * receiver, const PrimitiveID & sender  , walberla::mpi::RecvBuffer & buffer ) = 0;

  virtual void communicateLocalEdgeToFace( const Edge * sender, Face * receiver ) = 0;

};


} // namespace hhg
