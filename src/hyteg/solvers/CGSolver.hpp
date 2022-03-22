/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "core/Abort.h"
#include "core/timing/TimingTree.h"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

template < class OperatorType >
class CGSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   /// The algorithm is copied from the book: "Finite Elements and Fast Iterative Solvers"
   /// Therefore the variables are named like the ones in the book
   CGSolver(
       const std::shared_ptr< PrimitiveStorage >& storage,
       uint_t                                     minLevel,
       uint_t                                     maxLevel,
       uint_t                                     maxIter        = std::numeric_limits< uint_t >::max(),
       real_t                                     tolerance      = 1e-16,
       std::shared_ptr< Solver< OperatorType > >  preconditioner = std::make_shared< IdentityPreconditioner< OperatorType > >() )
   : p_( "p", storage, minLevel, maxLevel )
   , z_( "z", storage, minLevel, maxLevel )
   , ap_( "ap", storage, minLevel, maxLevel )
   , r_( "r", storage, minLevel, maxLevel )
   , preconditioner_( preconditioner )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , printInfo_( false )
   , tolerance_( tolerance )
   , restartFrequency_( std::numeric_limits< uint_t >::max() )
   , maxIter_( maxIter )
   , timingTree_( storage->getTimingTree() )
   {
      if ( !std::is_same< FunctionType, typename OperatorType::dstType >::value )
      {
         WALBERLA_ABORT( "CGSolver does not work for Operator with different src and dst FunctionTypes" );
      }
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      if ( maxIter_ == 0 )
         return;

      // if ( x.isDummy() || b.isDummy() )
      //    return;

      timingTree_->start( "CG Solver" );

      p_.copyBoundaryConditionFromFunction( x );
      z_.copyBoundaryConditionFromFunction( x );
      ap_.copyBoundaryConditionFromFunction( x );
      r_.copyBoundaryConditionFromFunction( x );

      real_t prsold = 0;
      init( A, x, b, level, prsold );
      real_t res_start = std::sqrt( r_.dotGlobal( r_, level, flag_ ) );

      if ( res_start < tolerance_ )
      {
         if ( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[CG] converged" );
         }
         timingTree_->stop( "CG Solver" );
         return;
      }
      real_t pAp, alpha, rsnew, sqrsnew, prsnew, beta;

      for ( size_t i = 0; i < maxIter_; ++i )
      {
         A.apply( p_, ap_, level, flag_, Replace );
         pAp = p_.dotGlobal( ap_, level, flag_ );

         alpha = prsold / pAp;
         x.add( {alpha}, {p_}, level, flag_ );
         r_.add( {-alpha}, {ap_}, level, flag_ );
         rsnew   = r_.dotGlobal( r_, level, flag_ );
         sqrsnew = std::sqrt( rsnew );

         if ( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[CG] residual: " << sqrsnew );
         }

         if ( sqrsnew < tolerance_ )
         {
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[CG] converged after " << i << " iterations" );
            }
            break;
         }

         preconditioner_->solve( A, z_, r_, level );
         prsnew = r_.dotGlobal( z_, level, flag_ );
         beta   = prsnew / prsold;

         p_.assign( {1.0, beta}, {z_, p_}, level, flag_ );
         prsold = prsnew;

         if ( i > 0 && i % restartFrequency_ == 0 )
         {
            init( A, x, b, level, prsold );
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[CG] restarted" );
            }
         }
      }
      timingTree_->stop( "CG Solver" );
   }

   /// \brief Compute tridiagonal matrix associated with underlying Lanzos process
   ///
   /// This method can be used to employ the Conjugate Gradient algorithm to setup
   /// the tridiagonal matrix T for a symmetric positive definite operator/matrix
   /// associated with underlying Lanzos process. This matrix e.g. provides information
   /// on the extremal eigenvalues of the operator.
   ///
   /// The matrix is defined when performing m iteration steps as:
   ///
   /// \f[
   /// T_m = \left(\begin{array}{cccccc}
   ///  \displaystyle \frac{1}{\alpha_0} & \displaystyle \frac{\sqrt{\beta_1}}{\alpha_0}
   ///       & & & & \\[2ex]
   ///  \displaystyle \frac{\sqrt{\beta_1}}{\alpha_0} &
   ///       \displaystyle \frac{1}{\alpha_1} + \frac{\beta_1}{\alpha_0} &
   ///       \displaystyle \frac{\sqrt{\beta_2}}{\alpha_1} & & & \\[2ex]
   ///  & \ddots & \ddots & \ddots & \\
   ///  & & \ddots & \ddots & \displaystyle\frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}} \\
   ///  & & & \displaystyle\frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}} &
   ///        \displaystyle\frac{1}{\alpha_{m-1}} + \frac{\beta_{m-2}}{\alpha_{m-1}}
   ///     \end{array}\right)
   /// \f]
   ///
   /// here \f$\alpha_k\f$ is the step length CG takes into the current search direction, while
   /// \f$\beta_k\f$ is the scalar used to update the search direction itself. For further details
   /// see e.g. the book "Iterative methods for sparse linear systems" by Yousef Saad.
   ///
   /// \param A          operator to be used in CG/Lanczos method
   /// \param x          auxilliary vector needed for performing CG iterations
   /// \param b          right-hand side vector used for CG iterations
   /// \param level      grid level to work on (operator & gridfunctions)
   /// \param numSteps   number of CG steps performed corresponds to dimension of matrix
   /// \param mainDiag   on return this vector containes the entries of T on the main diagonal
   /// \param subDiag    on return this vector containes the entries of T on the 1st sub-diagonal
   void setupLanczosTriDiagMatrix( const OperatorType&    A,
                                   const FunctionType&    x,
                                   const FunctionType&    b,
                                   const uint_t           level,
                                   const uint_t           numSteps,
                                   std::vector< real_t >& mainDiag,
                                   std::vector< real_t >& subDiag )
   {
      real_t prsold, pAp, prsnew, alpha, alpha_old, beta;

      // ----------------
      //  initialisation
      // ----------------

      // prepare vectors for Lanczos matrix data
      mainDiag.clear();
      mainDiag.reserve( numSteps );
      subDiag.clear();
      subDiag.reserve( numSteps - 1 );

      // init CG
      A.apply( x, p_, level, flag_, Replace );
      r_.assign( {1.0, -1.0}, {b, p_}, level, flag_ );
      p_.assign( {1.0}, {r_}, level, flag_ );
      prsold = r_.dotGlobal( r_, level, flag_ );

      // required for diagonal entries, set values
      // such that (1,1) entry is computed corretly
      alpha_old = 1.0;
      beta      = 0.0;

      // ---------------
      //  CG iterations
      // ---------------
      for ( uint_t i = 1; i < numSteps; ++i )
      {
         A.apply( p_, ap_, level, flag_, Replace );
         pAp = p_.dotGlobal( ap_, level, flag_ );

         alpha = prsold / pAp;
         mainDiag.push_back( 1.0 / alpha + beta / alpha_old );

         x.add( {alpha}, {p_}, level, flag_ );
         r_.add( {-alpha}, {ap_}, level, flag_ );

         prsnew = r_.dotGlobal( r_, level, flag_ );
         beta   = prsnew / prsold;
         subDiag.push_back( std::sqrt( beta ) / alpha );

         p_.assign( {1.0, beta}, {r_, p_}, level, flag_ );
         prsold = prsnew;

         alpha_old = alpha;
      }

      // final diagonal matrix entry
      A.apply( p_, ap_, level, flag_, Replace );
      pAp = p_.dotGlobal( ap_, level, flag_ );

      alpha = prsold / pAp;
      mainDiag.push_back( 1.0 / alpha + beta / alpha_old );
   }

   void setPrintInfo( bool printInfo ) { printInfo_ = printInfo; }

 private:
   void init( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level, real_t& prsold ) const
   {
      A.apply( x, p_, level, flag_, Replace );
      r_.assign( {1.0, -1.0}, {b, p_}, level, flag_ );
      preconditioner_->solve( A, z_, r_, level );
      p_.assign( {1.0}, {z_}, level, flag_ );
      prsold = r_.dotGlobal( z_, level, flag_ );
   }

   FunctionType                              p_;
   FunctionType                              z_;
   FunctionType                              ap_;
   FunctionType                              r_;
   std::shared_ptr< Solver< OperatorType > > preconditioner_;

   hyteg::DoFType flag_;
   bool           printInfo_;
   real_t         tolerance_;
   uint_t         restartFrequency_;
   uint_t         maxIter_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
};

} // namespace hyteg
