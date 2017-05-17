#ifndef FLAGS_HPP
#define FLAGS_HPP

namespace hhg
{

enum Boundary:size_t
{
  All = 1+2+4,
  Inner = 1,
  DirichletBoundary = 2,
  NeumannBoundary = 4,
};

inline Boundary operator|(Boundary a, Boundary b){
  return Boundary(a|b);
}
inline Boundary operator&(Boundary a, Boundary b){
  return Boundary(a&b);
}

inline bool testFlag(Boundary a, Boundary b)
{
  return (a & b) != 0;
}

}

#endif /* FLAGS_HPP */