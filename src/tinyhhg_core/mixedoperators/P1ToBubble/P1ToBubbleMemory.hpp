#pragma once

namespace hhg {

class FaceP1ToBubbleStencilMemory {
 public:

  typedef std::array<std::unique_ptr<real_t[]>, 2> StencilStack;
  std::map<size_t, StencilStack> data;

  inline StencilStack &addlevel(size_t level) {
    if (data.count(level) > 0)
    WALBERLA_LOG_WARNING("Level already exists.")
    else {
      data[level] = StencilStack{{ std::unique_ptr< real_t[ ]>( new real_t[ getGrayStencilSize(level) ] ),
                                   std::unique_ptr< real_t[ ]>( new real_t[ getBlueStencilSize(level) ] ) }};
    }
    return data[level];
  }

  inline size_t getGrayStencilSize(size_t level) {
    return 3;
  }

  inline size_t getBlueStencilSize(size_t level) {
    return 3;
  }

};

}
