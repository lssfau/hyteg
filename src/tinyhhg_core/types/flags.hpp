#pragma once

namespace hhg
{

enum UpdateType
{
  Replace = 0,
  Add = 1
};

enum DoFType:size_t
{
  All = 1+2+4+8,
  Boundary = 2+4+8,
  Inner = 1,
  DirichletBoundary = 2,
  NeumannBoundary = 4,
  FreeslipBoundary = 8
};

inline DoFType operator|(DoFType a, DoFType b){
  return DoFType(static_cast<size_t>(a) | static_cast<size_t>(b));
}

inline DoFType operator&(DoFType a, DoFType b){
  return DoFType(static_cast<size_t>(a) & static_cast<size_t>(b));
}

inline std::ostream& operator<<(std::ostream &os, const DoFType type){
  switch(type) {
    case FreeslipBoundary   :
      return os << "FreeslipBoundary";
    case NeumannBoundary   :
      return os << "NeumannBoundary";
    case DirichletBoundary :
      return os << "DirichletBoundary";
    default:
      return os << "Inner";
  }
}


inline bool testFlag(DoFType a, DoFType b)
{
  return (a & b) != 0;
}

enum class MemoryType {
  Base,
  P1Stencil,
  P1Function,
  P1BubbleFunction,
  P1BubbleStencil
};

} // namespace hhg

