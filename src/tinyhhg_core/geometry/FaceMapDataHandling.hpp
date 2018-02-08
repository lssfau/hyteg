#pragma once

#include "FaceMapMemory.hpp"

namespace hhg {

class MacroFaceFaceMapDataHandling : public OnlyInitializeDataHandling< MacroFaceFaceMap2DMemory, Face >
{
public:
  std::shared_ptr< MacroFaceFaceMap2DMemory > initialize( const Face * const face ) const {
    return std::make_shared<MacroFaceFaceMap2DMemory>();
  }
};

}