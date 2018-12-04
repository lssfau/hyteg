#pragma once

#include "core/DataTypes.h"

namespace hhg {
template < class OperatorType >
class Solver
{
 public:
   /// solves the system A * x = b
   virtual void solve( const OperatorType&             A,
                       const typename OperatorType::srcType& x,
                       const typename OperatorType::dstType& b,
                       const walberla::uint_t         level )const = 0;

   virtual ~Solver() = default;
};
} // namespace hhg