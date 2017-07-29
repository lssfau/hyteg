#pragma once

namespace hhg {

class VertexBubbleToP1StencilMemory {
 public:

  std::map<size_t, std::unique_ptr<real_t[]>> data;
  size_t num_deps_;

  inline std::unique_ptr<real_t[]> &addlevel(size_t level, size_t num_deps) {
    if (data.count(level) > 0)
      WALBERLA_LOG_WARNING("Level already exists.")
    else {
      this->num_deps_ = num_deps;
      data[level] = hhg::make_unique<real_t[]>(getSize(level));
    }
    return data[level];
  }

  inline size_t getSize(size_t level) {
    return num_deps_;
  }

};

class EdgeBubbleToP1StencilMemory {
 public:

  std::map<size_t, std::unique_ptr<real_t[]>> data;

  inline std::unique_ptr<real_t[]> &addlevel(size_t level) {
    if (data.count(level) > 0)
    WALBERLA_LOG_WARNING("Level already exists.")
    else {
      data[level] = hhg::make_unique<real_t[]>(getSize(level));
    }
    return data[level];
  }

  inline size_t getSize(size_t level) {
    return 6;
  }

};

class FaceBubbleToP1StencilMemory {
 public:

  std::map<size_t, std::unique_ptr<real_t[]>> data;

  inline std::unique_ptr<real_t[]> &addlevel(size_t level) {
    if (data.count(level) > 0)
    WALBERLA_LOG_WARNING("Level already exists.")
    else {
      data[level] = hhg::make_unique<real_t[]>(getSize(level));
    }
    return data[level];
  }

  inline size_t getSize(size_t level) {
    return 6;
  }

};

}