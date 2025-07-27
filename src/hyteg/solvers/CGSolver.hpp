/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Marcus Mohr, Nils Kohl, Andreas Burkhart.
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

#include "hyteg/functions/FunctionTools.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hyteg {

using walberla::uint_t;

template < class OperatorType >
class CGSolver : public Solver< OperatorType >
{
 public:
   using FunctionType = typename OperatorType::srcType;
   using ValueType    = typename FunctionTrait< FunctionType >::ValueType;

   /// The algorithm is copied from the book: "Finite Elements and Fast Iterative Solvers"
   /// Therefore the variables are named like the ones in the book
   CGSolver(
       const std::shared_ptr< PrimitiveStorage >& storage,
       uint_t                                     minLevel,
       uint_t                                     maxLevel,
       uint_t                                     maxIter           = std::numeric_limits< uint_t >::max(),
       typename FunctionType::valueType           relativeTolerance = 1e-16,
       typename FunctionType::valueType           absuluteTolerance = 1e-16,
       std::shared_ptr< Solver< OperatorType > >  preconditioner = std::make_shared< IdentityPreconditioner< OperatorType > >(),
       bool                                       lowMemoryMode  = false,
       std::function<
           bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t ) > iterationHook =
           []( uint_t it, const OperatorType& A, const FunctionType& x, const FunctionType& b, uint_t level, real_t initRes ) {
              WALBERLA_UNUSED( it );
              WALBERLA_UNUSED( A );
              WALBERLA_UNUSED( x );
              WALBERLA_UNUSED( b );
              WALBERLA_UNUSED( level );
              WALBERLA_UNUSED( initRes );
              return false;
           } )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , preconditioner_( preconditioner )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , printInfo_( false )
   , absuluteTolerance_( absuluteTolerance )
   , relativeTolerance_( relativeTolerance )
   , restartFrequency_( std::numeric_limits< uint_t >::max() )
   , maxIter_( maxIter )
   , iterations_( maxIter_ )
   , name_( "CG" )
   , timingTree_( storage->getTimingTree() )
   , iterationHook_( iterationHook )
   , lowMemoryMode_( lowMemoryMode )
   {
      if ( !std::is_same< FunctionType, typename OperatorType::dstType >::value )
      {
         WALBERLA_ABORT( "CGSolver does not work for Operator with different src and dst FunctionTypes" );
      }

      if ( !lowMemoryMode_ )
      {
         p_  = std::make_shared< FunctionType >( "p", storage, minLevel, maxLevel );
         z_  = std::make_shared< FunctionType >( "z", storage, minLevel, maxLevel );
         ap_ = std::make_shared< FunctionType >( "ap", storage, minLevel, maxLevel );
         r_  = std::make_shared< FunctionType >( "r", storage, minLevel, maxLevel );
      }
   }

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      if ( maxIter_ == 0 )
         return;

      timingTree_->start( "CG Solver" );

      std::shared_ptr< FunctionType > pSolve;
      std::shared_ptr< FunctionType > zSolve;
      std::shared_ptr< FunctionType > apSolve;
      std::shared_ptr< FunctionType > rSolve;

      if ( !lowMemoryMode_ )
      {
         pSolve  = p_;
         zSolve  = z_;
         apSolve = ap_;
         rSolve  = r_;
      }
      else
      {
         pSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         zSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         apSolve = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         rSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
      }

      copyBCs( x, *pSolve );
      copyBCs( x, *zSolve );
      copyBCs( x, *apSolve );
      copyBCs( x, *rSolve );

      pSolve->setToZero( level );
      zSolve->setToZero( level );
      apSolve->setToZero( level );
      rSolve->setToZero( level );

      typename FunctionType::valueType prsold = 0;
      init( A, x, b, level, prsold, pSolve, zSolve, rSolve );
      typename FunctionType::valueType res_start = std::sqrt( rSolve->dotGlobal( *rSolve, level, flag_ ) );

      if ( res_start < absuluteTolerance_ )
      {
         if ( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] converged" );
         }
         timingTree_->stop( "CG Solver" );
         return;
      }
      typename FunctionType::valueType pAp, alpha, rsnew, sqrsnew, prsnew, beta, relRes;

      for ( size_t i = 0; i < maxIter_; ++i )
      {
         A.apply( *pSolve, *apSolve, level, flag_, Replace );
         pAp = pSolve->dotGlobal( *apSolve, level, flag_ );

         alpha = prsold / pAp;
         x.add( { alpha }, { *pSolve }, level, flag_ );
         rSolve->add( { -alpha }, { *apSolve }, level, flag_ );
         rsnew   = rSolve->dotGlobal( *rSolve, level, flag_ );
         sqrsnew = std::sqrt( rsnew );
         relRes  = sqrsnew / res_start;

         if ( iterationHook_( i, A, x, b, level, res_start ) )
         {
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] iteration hook stopped the solver after " << std::defaultfloat << i
                                              << " iterations" );
            }
            iterations_ = i;
            break;
         }

         if ( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] iter: " << i << ", residual: " << sqrsnew
                                           << " ; relative residual: " << relRes );
         }

         if ( relRes < relativeTolerance_ || sqrsnew < absuluteTolerance_ )
         {
            iterations_ = i;
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] converged after " << i << " iterations" );
            }
            break;
         }

         zSolve->interpolate( walberla::numeric_cast< ValueType >( 0 ), level, All );
         preconditioner_->solve( A, *zSolve, *rSolve, level );
         prsnew = rSolve->dotGlobal( *zSolve, level, flag_ );
         beta   = prsnew / prsold;

         pSolve->assign( { walberla::numeric_cast< ValueType >( 1.0 ), beta }, { *zSolve, *pSolve }, level, flag_ );
         prsold = prsnew;

         if ( i > 0 && i % restartFrequency_ == 0 )
         {
            init( A, x, b, level, prsold, pSolve, zSolve, rSolve );
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[" << name_ << "] restarted" );
            }
         }
      }
      timingTree_->stop( "CG Solver" );
   }

   /** \brief Compute tridiagonal matrix associated with underlying Lanzos process
    *
    *  This method can be used to employ the Conjugate Gradient algorithm to setup
    *  the tridiagonal matrix T for a symmetric positive definite operator/matrix
    *  associated with underlying Lanzos process. This matrix e.g. provides information
    *  on the extremal eigenvalues of the operator.
    *
    *  The matrix is defined when performing m iteration steps as:
    *
    *  \f[
    *  T_m = \left(\begin{array}{cccccc}
    *   \displaystyle \frac{1}{\alpha_0} & \displaystyle \frac{\sqrt{\beta_1}}{\alpha_0}
    *        & & & & \\[2ex]
    *   \displaystyle \frac{\sqrt{\beta_1}}{\alpha_0} &
    *        \displaystyle \frac{1}{\alpha_1} + \frac{\beta_1}{\alpha_0} &
    *        \displaystyle \frac{\sqrt{\beta_2}}{\alpha_1} & & & \\[2ex]
    *   & \ddots & \ddots & \ddots & \\
    *   & & \ddots & \ddots & \displaystyle\frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}} \\
    *   & & & \displaystyle\frac{\sqrt{\beta_{m-2}}}{\alpha_{m-2}} &
    *         \displaystyle\frac{1}{\alpha_{m-1}} + \frac{\beta_{m-2}}{\alpha_{m-1}}
    *      \end{array}\right)
    *  \f]
    *
    *  here \f$\alpha_k\f$ is the step length CG takes into the current search direction, while
    *  \f$\beta_k\f$ is the scalar used to update the search direction itself. For further details
    *  see e.g. the book "Iterative methods for sparse linear systems" by Yousef Saad.
    *
    *  \param A          operator to be used in CG/Lanczos method
    *  \param x          auxilliary vector needed for performing CG iterations
    *  \param b          right-hand side vector used for CG iterations
    *  \param level      grid level to work on (operator & gridfunctions)
    *  \param numSteps   number of CG steps performed corresponds to dimension of matrix
    *  \param mainDiag   on return this vector containes the entries of T on the main diagonal
    *  \param subDiag    on return this vector containes the entries of T on the 1st sub-diagonal
   **/
   void setupLanczosTriDiagMatrix( const OperatorType&                              A,
                                   const FunctionType&                              x,
                                   const FunctionType&                              b,
                                   const uint_t                                     level,
                                   const uint_t                                     numSteps,
                                   std::vector< typename FunctionType::valueType >& mainDiag,
                                   std::vector< typename FunctionType::valueType >& subDiag )
   {
      std::shared_ptr< FunctionType > pSolve;
      std::shared_ptr< FunctionType > apSolve;
      std::shared_ptr< FunctionType > rSolve;

      if ( !lowMemoryMode_ )
      {
         pSolve  = p_;
         apSolve = ap_;
         rSolve  = r_;
      }
      else
      {
         pSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         apSolve = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );
         rSolve  = getTemporaryFunction< FunctionType >( storage_, minLevel_, maxLevel_ );

         pSolve->setToZero( level );
         apSolve->setToZero( level );
         rSolve->setToZero( level );
      }

      typename FunctionType::valueType prsold, pAp, prsnew, alpha, alpha_old, beta;

      // ----------------
      //  initialisation
      // ----------------

      // prepare vectors for Lanczos matrix data
      mainDiag.clear();
      mainDiag.reserve( numSteps );
      subDiag.clear();
      subDiag.reserve( numSteps - 1 );

      // init CG
      A.apply( x, *pSolve, level, flag_, Replace );
      rSolve->assign( { walberla::numeric_cast< ValueType >( 1.0 ), walberla::numeric_cast< ValueType >( -1.0 ) },
                      { b, *pSolve },
                      level,
                      flag_ );
      pSolve->assign( { walberla::numeric_cast< ValueType >( 1.0 ) }, { *rSolve }, level, flag_ );
      prsold = rSolve->dotGlobal( *rSolve, level, flag_ );

      // required for diagonal entries, set values
      // such that (1,1) entry is computed corretly
      alpha_old = 1.0;
      beta      = 0.0;

      // ---------------
      //  CG iterations
      // ---------------
      for ( uint_t i = 1; i < numSteps; ++i )
      {
         A.apply( *pSolve, *apSolve, level, flag_, Replace );
         pAp = pSolve->dotGlobal( *apSolve, level, flag_ );

         alpha = prsold / pAp;
         mainDiag.push_back( 1.0 / alpha + beta / alpha_old );

         x.add( { alpha }, { *pSolve }, level, flag_ );
         rSolve->add( { -alpha }, { *apSolve }, level, flag_ );

         prsnew = rSolve->dotGlobal( *rSolve, level, flag_ );
         beta   = prsnew / prsold;
         subDiag.push_back( std::sqrt( beta ) / alpha );

         pSolve->assign( { walberla::numeric_cast< ValueType >( 1.0 ), beta }, { *rSolve, *pSolve }, level, flag_ );
         prsold = prsnew;

         alpha_old = alpha;
      }

      // final diagonal matrix entry
      A.apply( *pSolve, *apSolve, level, flag_, Replace );
      pAp = pSolve->dotGlobal( *apSolve, level, flag_ );

      alpha = prsold / pAp;
      mainDiag.push_back( real_c( 1.0 ) / alpha + beta / alpha_old );
   }

   uint_t getIterations() { return iterations_; }

   void setPrintInfo( bool printInfo ) { printInfo_ = printInfo; }
   void setName( std::string newName ) { name_ = newName; }
   void setDoFType( hyteg::DoFType flag ) { flag_ = flag; }
   void setIterationHook(
       std::function< bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t ) >
           iterationHook )
   {
      iterationHook_ = iterationHook;
   };

 private:
   void init( const OperatorType&                    A,
              const FunctionType&                    x,
              const FunctionType&                    b,
              const uint_t                           level,
              typename FunctionType::valueType&      prsold,
              const std::shared_ptr< FunctionType >& pSolve,
              const std::shared_ptr< FunctionType >& zSolve,
              const std::shared_ptr< FunctionType >& rSolve ) const
   {
      A.apply( x, *pSolve, level, flag_, Replace );
      rSolve->assign( { walberla::numeric_cast< ValueType >( 1.0 ), walberla::numeric_cast< ValueType >( -1.0 ) },
                      { b, *pSolve },
                      level,
                      flag_ );
      zSolve->interpolate( walberla::numeric_cast< ValueType >( 0 ), level, All );
      preconditioner_->solve( A, *zSolve, *rSolve, level );
      pSolve->assign( { walberla::numeric_cast< ValueType >( 1.0 ) }, { *zSolve }, level, flag_ );
      prsold = rSolve->dotGlobal( *zSolve, level, flag_ );
   }
   std::shared_ptr< hyteg::PrimitiveStorage > storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;

   std::shared_ptr< FunctionType >           p_;
   std::shared_ptr< FunctionType >           z_;
   std::shared_ptr< FunctionType >           ap_;
   std::shared_ptr< FunctionType >           r_;
   std::shared_ptr< Solver< OperatorType > > preconditioner_;

   hyteg::DoFType                   flag_;
   bool                             printInfo_;
   typename FunctionType::valueType absuluteTolerance_;
   typename FunctionType::valueType relativeTolerance_;
   uint_t                           restartFrequency_;
   uint_t                           maxIter_;
   uint_t                           iterations_;

   std::string name_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   std::function< bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t ) > iterationHook_;

   bool lowMemoryMode_;
};

} // namespace hyteg
