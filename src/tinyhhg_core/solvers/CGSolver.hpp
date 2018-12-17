#pragma once

#include "core/Abort.h"

#include "tinyhhg_core/solvers/Solver.hpp"
#include "tinyhhg_core/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hhg {

using walberla::real_t;
using walberla::uint_t;

template < class OperatorType >
class CGSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

//   static_assert( !std::is_same< FunctionType, typename OperatorType::dstType >::value,
//                  "CGSolver does not work for Operator with different src and dst FunctionTypes" );

   /// The algorithm is copied from the book: "Finite Elements and Fast Iterative Solvers"
   /// Therefore the variables are named like the ones in the book
   CGSolver( const std::shared_ptr< PrimitiveStorage >& storage,
             uint_t                                     minLevel,
             uint_t                                     maxLevel,
             uint_t                                     maxIter = std::numeric_limits< uint_t>::max(),
             real_t                                     tolerance = 1e-16,
             std::shared_ptr< Solver< OperatorType > > preconditioner = std::make_shared< IdentityPreconditioner< OperatorType > >() )
   : p_( "p", storage, minLevel, maxLevel )
   , z_( "z", storage, minLevel, maxLevel )
   , ap_( "ap", storage, minLevel, maxLevel )
   , r_( "r", storage, minLevel, maxLevel )
   , preconditioner_( preconditioner )
   , flag_( hhg::Inner | hhg::NeumannBoundary )
   , printInfo_( false )
   , tolerance_( tolerance )
   , restartFrequency_( std::numeric_limits< uint_t >::max() )
   , maxIter_( maxIter )
   {
      if( !std::is_same< FunctionType, typename OperatorType::dstType >::value )
      {
         WALBERLA_ABORT( "CGSolver does not work for Operator with different src and dst FunctionTypes" );
      }
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      real_t prsold = 0;
      init( A, x, b, level, prsold );
      real_t res_start = std::sqrt( r_.dotGlobal( r_, level, flag_ ) );

      if( res_start < tolerance_ )
      {
         if( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[CG] converged" );
         }
         return;
      }
      real_t pAp, alpha, rsnew, sqrsnew, prsnew, beta;

      for( size_t i = 0; i < maxIter_; ++i )
      {
         A.apply( p_, ap_, level, flag_, Replace );
         pAp = p_.dotGlobal( ap_, level, flag_ );

         alpha = prsold / pAp;
         x.add( {alpha}, {p_}, level, flag_ );
         r_.add( {-alpha}, {ap_}, level, flag_ );
         rsnew   = r_.dotGlobal( r_, level, flag_ );
         sqrsnew = std::sqrt( rsnew );

         if( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[CG] residual: " << sqrsnew );
         }

         if( sqrsnew / res_start < tolerance_ )
         {
            if( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[CG] converged after " << i << " iterations" );
            }
            break;
         }

         preconditioner_->solve( A, r_, z_, level );
         prsnew = r_.dotGlobal( z_, level, flag_ );
         beta   = prsnew / prsold;

         p_.assign( {1.0, beta}, {z_, p_}, level, flag_ );
         prsold = prsnew;

         if( i % restartFrequency_ == 0 )
         {
            init( A, x, b, level, prsold );
            if( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[CG] restarted" );
            }
         }
      }
   }

 private:
   void init( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level, real_t& prsold ) const
   {
      A.apply( x, p_, level, flag_, Replace );
      r_.assign( {1.0, -1.0}, {b, p_}, level, flag_ );
      preconditioner_->solve( A, r_, z_, level );
      p_.assign( {1.0}, {z_}, level, flag_ );
      prsold = r_.dotGlobal( z_, level, flag_ );
   }

   FunctionType                          p_;
   FunctionType                          z_;
   FunctionType                          ap_;
   FunctionType                          r_;
   std::shared_ptr< Solver< OperatorType > > preconditioner_;

   hhg::DoFType flag_;
   bool         printInfo_;
   real_t       tolerance_;
   uint_t       restartFrequency_;
   uint_t       maxIter_;
};

} // namespace hhg
