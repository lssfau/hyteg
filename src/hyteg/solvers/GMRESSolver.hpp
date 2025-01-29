/*
* Copyright (c) 2024 Maximilian Dechant, Andreas Wagner, Andreas Burkhart.
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

#include <fstream>
#include <iostream>

#include "core/Abort.h"
#include "core/timing/TimingTree.h"

#include "hyteg/functions/FunctionTools.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

#include "../eigen/Eigen/Eigen"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// Left preconditioned FGMRES implementation following:
/// Saad, Youcef and Schultz, Martin H.
/// "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems "
/// SIAM Journal on Scientific and Statistical Computing , Vol. 7, No. 3, p. 856-869, 1986
/// https://doi.org/10.1137/0907058
///
/// \param maxKrylowDim Maximum number of iterations, i.e. maximum krylov dimension
/// \param restartLength Restart the algorithm after restartLength steps. Restarting means, that we discard all data from previous steps (especially the saved functions) and start from the previous approximate solution as an initial guess.
/// \param arnoldiTOL The loop constructing an orthogonal basis for the preconditioned Krylov space is called "Arnoldi loop" or "Arnoldi process". If the 2-norm of a newly formed basis vector (before normalisation) is smaller than the arnoldi tolerance, the algorithm terminates.
/// \param approxTOL Termination tolerance for the approximated residual
/// \param doubleOrthoTOL Due to rounding errors it is sometimes necessary to repeat the orthogonalisation process. If doubleOrthoTOL == 0 or if the 2-Norm of the difference between the orthogonalised and non-orthogonalised vector is greater than doubleOrthoTOL (calculation required) a second orthogonalisation process is triggered. Since the double orthogonalisation step is most of the cheaper than computing the orthogonalisation difference (preconditioner solve required) you should leave this parameter at 0 if you do not have a specific reason to do otherwise.
/// \param preconditioner A given left preconditioner
/// \param iterationHook An iteration hook that is called after every iteration. The iteration hook is expected to return true if you want to stop the solver loop. The hook can for example be used to calculate the L2 error after each iteration.
/// \param generateSolutionForHook Defines wether the iteration hook is provided the current approximate solution calculated from the saved basis functions. Otherwise the x input of the iteration hook should not be used. Calculating the approximate solution from the given basis functions is costly, so only use this if necessary.
template < class OperatorType >
class GMRESSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   GMRESSolver(
       const std::shared_ptr< PrimitiveStorage >& storage,
       uint_t                                     minLevel,
       uint_t                                     maxLevel,
       uint_t                                     maxKrylowDim   = 1000,
       uint_t                                     restartLength  = 1000,
       real_t                                     arnoldiTOL     = 1e-16,
       real_t                                     approxTOL      = 1e-16,
       real_t                                     doubleOrthoTOL = 0,
       std::shared_ptr< Solver< OperatorType > >  preconditioner = std::make_shared< IdentityPreconditioner< OperatorType > >(),
       std::function< bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t, real_t ) >
           iterationHook =
               []( uint_t              it,
                   const OperatorType& A,
                   const FunctionType& x,
                   const FunctionType& b,
                   uint_t              level,
                   real_t              ApproxError,
                   real_t              wNorm ) {
                  WALBERLA_UNUSED( it );
                  WALBERLA_UNUSED( A );
                  WALBERLA_UNUSED( x );
                  WALBERLA_UNUSED( b );
                  WALBERLA_UNUSED( level );
                  WALBERLA_UNUSED( ApproxError );
                  WALBERLA_UNUSED( wNorm );
                  return false;
               },
       bool generateSolutionForHook = false )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , maxKrylowDim_( maxKrylowDim )
   , restartLength_( restartLength )
   , numberOfIterations_( maxKrylowDim )
   , arnoldiTOL_( arnoldiTOL )
   , approxTOL_( approxTOL )
   , doubleOrthoTOL_( doubleOrthoTOL )
   , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , printInfo_( false )
   , preconditioner_( preconditioner )
   , r0_( "r0", storage_, minLevel_, maxLevel_ )
   , rPrec_( "rPrec", storage_, minLevel_, maxLevel_ )
   , wPrec_( "wPrec", storage_, minLevel_, maxLevel_ )
   , orthoDiff_( "orthoDiff", storage_, minLevel_, maxLevel_ )
   , timingTree_( storage->getTimingTree() )
   , iterationHook_( iterationHook )
   , generateSolutionForHook_( generateSolutionForHook )
   {}

   uint_t getIterations() { return numberOfIterations_; }

   void setPrintInfo( bool printInfo ) { printInfo_ = printInfo; }

   void setMaxIter( uint_t maxIter ) { maxKrylowDim_ = maxIter; }

   void setRestartLength( uint_t restartLength ) { restartLength_ = restartLength; }

   void setPreconditioner( std::shared_ptr< Solver< OperatorType > > preconditioner ) { preconditioner_ = preconditioner; }

   void setAbsoluteTolerance( double atol ) { approxTOL_ = atol; }

   void setDoubleOrthogonalizationTolerance( double tol ) { doubleOrthoTOL_ = tol; }

   ~GMRESSolver() = default;

   void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
   {
      timingTree_->start( "GMRES Solver" );

      copyBCs( x, r0_ );
      copyBCs( x, rPrec_ );
      copyBCs( x, wPrec_ );
      copyBCs( x, orthoDiff_ );

      real_t approxERR = approxTOL_ + 1;
      bool   callback  = false;

      init( A, x, b, level );

      for ( uint_t j = 1; j < maxKrylowDim_; j++ )
      {
         // handle restart if applicable
         if ( j % restartLength_ == 0 )
         {
            generateFinalApproximation( x, level );
            init( A, x, b, level );
            WALBERLA_LOG_INFO_ON_ROOT( "[GMRES] restarted " );
            continue;
         }
         int currentIndex = j % restartLength_;

         // vector storage housekeeping
         if ( vecV_.size() <= currentIndex )
         {
            FunctionType w( "w", storage_, minLevel_, maxLevel_ );
            copyBCs( x, w );
            vecV_.push_back( w );
         }
         else
         {
            vecV_[currentIndex].interpolate( 0.0, level, flag_ );
         }

         // main algorithmic steps (b), 1 and (b), 2
         A.apply( vecV_[currentIndex - 1], vecV_[currentIndex], level, flag_ ); // (b), 1, left preconditioned
         wPrec_.interpolate( 0.0, level, flag_ );                               // (b), 1, left preconditioned
         preconditioner_->solve( A, wPrec_, vecV_[currentIndex], level );       // (b), 1, left preconditioned
         vecV_[currentIndex].assign( { 1.0 }, { wPrec_ }, level, flag_ );       // (b), 1, left preconditioned

         // resize H to the appropriate size
         H_.conservativeResize( currentIndex + 1, currentIndex );
         H_.block( currentIndex, 0, 1, currentIndex ) = MatrixXr::Zero(1, currentIndex);
         H_.block( 0, currentIndex-1, currentIndex + 1, 1 ) = MatrixXr::Zero(currentIndex + 1, 1);

         // (b), 3
         for ( uint_t i = 1; i <= currentIndex; i++ )
         {
            real_t h                      = vecV_[currentIndex].dotGlobal( vecV_[i - 1], level, flag_ );
            H_( i - 1, currentIndex - 1 ) = h;
            vecV_[currentIndex].add( { -h }, { vecV_[i - 1] }, level, flag_ );
         }

         // check if double orthogonalisation should be used
         if ( doubleOrthoTOL_ > 0 )
         {
            A.apply( vecV_[currentIndex - 1], orthoDiff_, level, flag_ );
            wPrec_.interpolate( 0.0, level, flag_ );
            preconditioner_->solve( A, wPrec_, orthoDiff_, level );
            orthoDiff_.assign( { 1.0, -1.0 }, { wPrec_, vecV_[currentIndex] }, level, flag_ );
         }

         if ( doubleOrthoTOL_ <= 0 || std::sqrt( orthoDiff_.dotGlobal( orthoDiff_, level, flag_ ) ) > doubleOrthoTOL_ )
         {
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[GMRES] invoked double-orthogonalization at iteration " << j );
            }
            for ( uint_t i = 1; i <= ( currentIndex ); i++ )
            {
               real_t h = vecV_[currentIndex].dotGlobal( vecV_[i - 1], level, flag_ );
               H_( i - 1, currentIndex - 1 ) += h;
               vecV_[currentIndex].add( { -h }, { vecV_[i - 1] }, level, flag_ );
            }
         }

         real_t wNorm = std::sqrt( vecV_[currentIndex].dotGlobal( vecV_[currentIndex], level, flag_ ) ); // (b), 4
         H_( currentIndex, currentIndex - 1 ) = wNorm;                                                   // (b), 4
         y_                                   = hessenbergMinimizer( beta_, H_, Q_, 1, approxERR );      // compute y_m
         if ( printInfo_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "[GMRES] approximated residual after " << j << " iterations : " << approxERR );
         }
         vecV_[currentIndex].assign( { real_c( 1.0 ) / wNorm }, { vecV_[currentIndex] }, level, flag_ ); // (b), 4

         if ( generateSolutionForHook_ )
         {
            // reuse r0 to hold the current solution
            r0_.assign( { 1.0 }, { x }, level, All );
            generateFinalApproximation( r0_, level );
            callback = iterationHook_( j, A, r0_, b, level, approxERR, wNorm );
         }
         else
         {
            callback = iterationHook_( j, A, x, b, level, approxERR, wNorm );
         }

         if ( callback )
         {
            if ( printInfo_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "[GMRES] iteration hook stopped the solver after " << std::defaultfloat << j
                                                                                             << " iterations" );
            }
            numberOfIterations_ = j;
            break;
         }

         if ( wNorm <= arnoldiTOL_ || approxERR <= approxTOL_ )
         {
            numberOfIterations_ = j;
            break;
         }
      }
      generateFinalApproximation( x, level );

      timingTree_->stop( "GMRES Solver" );
      return;
   }

 private:
   void init( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level )
   {
      // start of the algorithm
      A.apply( x, r0_, level, flag_ );                       // (a), 1
      r0_.assign( { 1.0, -1.0 }, { b, r0_ }, level, flag_ ); // (a), 1

      // compute the norm of the initial residual
      rPrec_.interpolate( 0.0, level, flag_ );                 // (a), 2, left preconditioned residual
      preconditioner_->solve( A, rPrec_, r0_, level );         // (a), 2, left preconditioned residual
      r0_.assign( { 1.0 }, { rPrec_ }, level, flag_ );         // (a), 2, left preconditioned residual
      beta_ = std::sqrt( r0_.dotGlobal( r0_, level, flag_ ) ); // (a), 2, left preconditioned residual

      // init matrices H, Q
      H_ = MatrixXr::Zero( 0, 0 );
      Q_ = MatrixXr::Ones( 1, 1 );

      // vector storage housekeeping
      if ( vecV_.empty() )
      {
         FunctionType v0( "v0", storage_, minLevel_, maxLevel_ );
         copyBCs( x, v0 );
         v0.assign( { real_c( 1.0 ) / beta_ }, { r0_ }, level, flag_ );
         vecV_.push_back( v0 );
      }
      else
      {
         vecV_[0].assign( { real_c( 1.0 ) / beta_ }, { r0_ }, level, flag_ );
      }
   }

   const std::shared_ptr< PrimitiveStorage >& storage_;
   uint_t                                     minLevel_;
   uint_t                                     maxLevel_;
   uint_t                                     maxKrylowDim_;
   uint_t                                     restartLength_;
   uint_t                                     numberOfIterations_;
   real_t                                     arnoldiTOL_;
   real_t                                     approxTOL_;
   real_t                                     doubleOrthoTOL_;
   hyteg::DoFType                             flag_;
   bool                                       printInfo_;
   std::shared_ptr< Solver< OperatorType > >  preconditioner_;

   real_t                      beta_;
   std::vector< FunctionType > vecV_;
   MatrixXr                    H_;
   MatrixXr                    Q_;
   VectorXr                    y_;
   FunctionType                r0_;
   FunctionType                rPrec_;
   FunctionType                wPrec_;
   FunctionType                orthoDiff_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;
   std::function< bool( uint_t, const OperatorType&, const FunctionType&, const FunctionType&, uint_t, real_t, real_t ) >
        iterationHook_;
   bool generateSolutionForHook_;

   void generateFinalApproximation( const FunctionType& x, uint_t level )
   {
      for ( int i = 0; i < y_.size(); i++ )
      {
         x.add( { y_( i ) }, { vecV_[i] }, level, flag_ );
      }
   }

   inline VectorXr getUnitVector( int length, int j, real_t val )
   {
      VectorXr answer = VectorXr::Zero( length );
      answer( j )     = val;
      return answer;
   }

   inline void expandQMatrix( MatrixXr& inputMatrix, int expandBy )
   {
      Eigen::Index originalRows = inputMatrix.rows();
      Eigen::Index originalCols = inputMatrix.cols();

      inputMatrix.conservativeResize( originalRows + expandBy, originalCols + expandBy );

      // new values need to be initialised
      inputMatrix.block( originalRows, 0, expandBy, originalCols ) = MatrixXr::Zero(expandBy, originalCols);
      inputMatrix.block( 0, originalCols, originalRows, expandBy ) = MatrixXr::Zero(originalRows, expandBy);
      inputMatrix.block( originalRows, originalCols, expandBy, expandBy ) = MatrixXr::Identity( expandBy, expandBy );

      // // Alternative:

      // MatrixXr outputMatrix = MatrixXr::Identity(inputMatrix.rows() + expandBy, inputMatrix.cols() + expandBy);

      // outputMatrix.block(0,0,inputMatrix.rows(),inputMatrix.cols()) = inputMatrix;

      // return outputMatrix;
   }

   VectorXr triangSolver( MatrixXr triMat, VectorXr targetVector )
   {
      VectorXr answer = VectorXr::Zero( triMat.cols() );
      for ( int i = triMat.cols() - 1; i >= 0; i-- )
      {
         real_t targetOffset = 0;
         for ( int j = triMat.cols() - 1; j > i; j-- )
         {
            targetOffset += triMat( i, j ) * answer( j );
         }
         answer( i ) = ( targetVector( i ) - targetOffset ) / triMat( i, i );
      }
      return answer;
   }

   VectorXr hessenbergMinimizer( real_t beta, MatrixXr& H, MatrixXr& Q, int numUnfinishedColumns, real_t& approxERR )
   {
      VectorXr approxVector;

      if ( H.rows() != Q.rows() )
      {
         expandQMatrix( Q, H.rows() - Q.rows() );
      }

      for ( int i = H.cols() - numUnfinishedColumns; i < H.cols(); i++ )
      {
         H.col( i ) = Q * H.col( i );
      }

      for ( int i = H.cols() - numUnfinishedColumns; i < H.cols(); i++ )
      {
         if ( H( i + 1, i ) == 0 )
         {
            continue;
         }
         real_t   s           = H( i + 1, i ) / std::sqrt( std::pow( H( i, i ), 2 ) + std::pow( H( i + 1, i ), 2 ) );
         real_t   c           = H( i, i ) / std::sqrt( std::pow( H( i, i ), 2 ) + std::pow( H( i + 1, i ), 2 ) );
         MatrixXr QNew        = MatrixXr::Identity( H.rows(), H.rows() );
         QNew( i, i )         = c;
         QNew( i + 1, i + 1 ) = c;
         QNew( i + 1, i )     = -s;
         QNew( i, i + 1 )     = s;

         H = QNew * H;
         Q = QNew * Q;
      }

      MatrixXr equationLeftSide  = H.block( 0, 0, H.cols(), H.cols() );
      VectorXr targetVector      = Q * getUnitVector( H.rows(), 0, beta );
      VectorXr equationRightSide = targetVector.head( H.cols() );
      approxERR                  = std::abs( targetVector( targetVector.rows() - 1 ) );

      return triangSolver( equationLeftSide, equationRightSide );
   }
};

} // namespace hyteg
