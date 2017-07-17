
#pragma once

namespace hhg {

using walberla::uint_t;

class NeighborPrimitive
{
public:

  uint_t getRank() const { return rank_; }

protected:

  NeighborPrimitive( const SetupPrimitive & setupPrimitive ) :
    rank_( setupPrimitive.getTargetRank() )
  {}

private:

  uint_t rank_;

};

}
