#pragma once

namespace hhg {

template < class OperatorType >
class GaussSeidelSmoother : public Solver< OperatorType >
{
 public:
   GaussSeidelSmoother()
   : flag_( hhg::Inner | hhg::NeumannBoundary )
   {}

   void solve( const OperatorType&             A,
               typename OperatorType::srcType& x,
               typename OperatorType::dstType& b,
               const uint_t&                   level ) override
   {
      A.smooth_gs( x, b, level, flag_ );
   }

 private:
   DoFType flag_;
};

} // namespace hhg