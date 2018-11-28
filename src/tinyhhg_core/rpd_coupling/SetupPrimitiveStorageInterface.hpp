
#pragma once

#include "rpd/domain/IDomain.h"

#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hhg {
namespace rpd {

using walberla::real_t;
using walberla::uint_t;
using walberla::rpd::Vec3;

class SetupPrimitiveStorageInterface : public walberla::rpd::IDomain
{
public:

  SetupPrimitiveStorageInterface( const std::shared_ptr< SetupPrimitiveStorage > & setupStorage ) :
    setupStorage_( setupStorage )
  {}

  virtual ~SetupPrimitiveStorageInterface() {}

  bool isContainedInProcessSubdomain( const uint_t rank, const Vec3& pt ) const override;

  int findContainingProcessRank( const Vec3& pt ) const override;

  void periodicallyMapToDomain( Vec3 & pt ) const override;

  std::vector<uint_t> getNeighborProcesses() const override;

  bool intersectsWithProcessSubdomain( const uint_t rank, const Vec3& pt, const real_t& radius ) const override;

  void correctParticlePosition( Vec3 & pt ) const override;

private:

  std::shared_ptr< SetupPrimitiveStorage > setupStorage_;
  std::vector< uint_t > neighborProcesses_;

};

}
}