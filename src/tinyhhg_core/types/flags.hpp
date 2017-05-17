#ifndef FLAGS_HPP
#define FLAGS_HPP

namespace hhg
{

enum UpdateType
{
  Replace = 0,
  Add = 1
};

enum DoFType:size_t
{
  All = 1+2+4,
  Inner = 1,
  DirichletBoundary = 2,
  NeumannBoundary = 4,
};

inline DoFType operator|(DoFType a, DoFType b){
  return DoFType(static_cast<size_t>(a) | static_cast<size_t>(b));
}

inline DoFType operator&(DoFType a, DoFType b){
  return DoFType(static_cast<size_t>(a) & static_cast<size_t>(b));
}

inline bool testFlag(DoFType a, DoFType b)
{
  return (a & b) != 0;
}

}

#endif /* FLAGS_HPP */
