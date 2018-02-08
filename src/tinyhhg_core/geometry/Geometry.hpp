#pragma once

#include "FaceMapDataHandling.hpp"

namespace hhg
{

class Geometry {
public:

  Geometry(const std::shared_ptr <PrimitiveStorage> &storage)
      : storage_(storage) {

    auto faceFaceMapMatrixMemoryDataHandling = std::make_shared< MacroFaceFaceMapDataHandling >();
    storage_->addFaceData(faceMapID_, faceFaceMapMatrixMemoryDataHandling, "MacroFaceFaceMapData");
  }

  const PrimitiveDataID<MacroFaceFaceMap2DMemory, Face>& getFaceMapID() const {
    return faceMapID_;
  };

private:
  std::shared_ptr <PrimitiveStorage> storage_;
  PrimitiveDataID<MacroFaceFaceMap2DMemory, Face> faceMapID_;
};

}