#ifndef FLAGS_HPP
#define FLAGS_HPP

namespace hhg
{

enum Boundary
{
  All = 1+2+4,
  Inner = 1,
  DirichletBoundary = 2,
  NeumannBoundary = 4,
};

inline bool testFlag(Boundary a, Boundary b)
{
  return (a & b);
}

}

#endif /* FLAGS_HPP */