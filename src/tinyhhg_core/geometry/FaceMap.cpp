#include "IdentityMap.hpp"
#include "CircularMap.hpp"
#include "AffineMap.hpp"

namespace hhg {

std::shared_ptr<FaceMap> FaceMap::deserialize(walberla::mpi::RecvBuffer& recvBuffer) {
  Type type;
  recvBuffer >> type;

  switch (type) {
    case IDENTITY:
      return std::shared_ptr<FaceMap>(new IdentityMap());
    case AFFINE:
      return std::shared_ptr<FaceMap>(new AffineMap(recvBuffer));
    case CIRCULAR:
      return std::shared_ptr<FaceMap>(new CircularMap(recvBuffer));
  }
}

}