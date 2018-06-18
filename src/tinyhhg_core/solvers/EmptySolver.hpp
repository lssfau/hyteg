
#pragma once

namespace hhg {

template< class Function_T, class Operator_T >
class EmptySolver
{
public:
    void solve(Operator_T &, Function_T &, Function_T &, Function_T &, uint_t, real_t, size_t, DoFType, bool)
    {
      // does nothing
    }
};

}