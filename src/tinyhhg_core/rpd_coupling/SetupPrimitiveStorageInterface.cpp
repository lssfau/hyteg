
#include "rpd/domain/IDomain.h"
#include "rpd/data/DataTypes.h"

#include "tinyhhg_core/rpd_coupling/SetupPrimitiveStorageInterface.hpp"
#include "tinyhhg_core/geometry/Intersection.hpp"

#include <algorithm>

namespace hhg {
namespace rpd {

using walberla::real_t;
using walberla::uint_t;
using walberla::rpd::Vec3;


bool SetupPrimitiveStorageInterface::isContainedInProcessSubdomain(const uint_t rank, const Vec3& pt) const
{
  const int containingRank = findContainingProcessRank( pt );
  return (int)rank == containingRank;
}

int SetupPrimitiveStorageInterface::findContainingProcessRank(const Vec3& pt) const
{
  std::set< int > containingProcessRanks;
  for ( const auto & cellIt : setupStorage_->getCells() )
  {
    auto cellID = cellIt.first;
    auto cell = cellIt.second;

    Point3D pointOfInterest({ pt[0], pt[1], pt[2] });

    if ( isPointInTetrahedron( pointOfInterest, cell->getCoordinates().at(0), cell->getCoordinates().at(1), cell->getCoordinates().at(2), cell->getCoordinates().at(3) ) )
    {
      containingProcessRanks.insert( static_cast< int >( setupStorage_->getTargetRank( cellID ) ) );
    }

  }
  if ( containingProcessRanks.empty() )
  {
    return -1;
  }
  else
  {
    // std::set is required to be sorted by the standard.
    // ::begin() returns the smallest element.
    return *containingProcessRanks.begin();
  }
}

void SetupPrimitiveStorageInterface::periodicallyMapToDomain(Vec3& pt) const
{
  WALBERLA_UNUSED( pt );
}

std::vector<uint_t> SetupPrimitiveStorageInterface::getNeighborProcesses() const
{
  std::set< uint_t > neighborProcesses;
  for ( const auto & cellIt : setupStorage_->getCells() )
  {
    auto cell = cellIt.second;

    for ( const auto & neighborVertexID : cell->neighborVertices() )
    {
      auto neighborVertex = setupStorage_->getVertex( neighborVertexID );
      for ( const auto neighborCellID : neighborVertex->neighborCells() )
      {
        neighborProcesses.insert( setupStorage_->getTargetRank( neighborCellID ) );
      }
    }
  }
  neighborProcesses.erase( uint_c( walberla::mpi::MPIManager::instance()->rank() ) );

  return std::vector< uint_t >( neighborProcesses.begin(), neighborProcesses.end() );
}

bool SetupPrimitiveStorageInterface::intersectsWithProcessSubdomain(const uint_t rank, const Vec3& pt, const real_t& radius) const
{
  for ( const auto & cellIt : setupStorage_->getCells() )
  {
    auto cellID = cellIt.first;
    auto cell = cellIt.second;

    Point3D pointOfInterest({ pt[0], pt[1], pt[2] });

    if ( setupStorage_->getTargetRank( cellID ) == rank
         && sphereTetrahedronIntersection( pointOfInterest, radius, cell->getCoordinates().at(0), cell->getCoordinates().at(1), cell->getCoordinates().at(2), cell->getCoordinates().at(3) ) )
    {
      return true;
    }

  }
  return false;
}

void SetupPrimitiveStorageInterface::correctParticlePosition(Vec3& pt) const
{
  WALBERLA_UNUSED( pt );
}

}
}
