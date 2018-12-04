#pragma once

namespace hhg {

template < class OperatorType >
class JacobiPreconditioner : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;
   JacobiPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, uint_t iterations )
   : iterations_( iterations )
   , tmp_( "jac_tmp", storage, minLevel, maxLevel )
   , flag_( hhg::Inner | hhg::NeumannBoundary )
   {}

   // y = M^{-1} * x
   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) const override
   {
      b.assign( {1.0}, {x}, level, flag_ );

      for( uint_t i = 0; i < iterations_; ++i )
      {
         tmp_.assign( {1.0}, {b}, level, flag_ );
         A.smooth_jac( b, x, tmp_, level, flag_ );
      }
   }

 private:
   uint_t       iterations_;
   FunctionType tmp_;
   hhg::DoFType flag_;
};

} // namespace hhg