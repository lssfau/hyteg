/*
 * Copyright (c) 2021 Maximilian Dechant.
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

#include "../eigen/Eigen/Eigen"

#include "core/Abort.h"
#include "core/timing/TimingTree.h"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

template < class OperatorType >
class GMRESSolver : public Solver< OperatorType >
{
 public:
   typedef typename OperatorType::srcType FunctionType;

   GMRESSolver(
      const std::shared_ptr< PrimitiveStorage >& storage,
      uint_t                                     minLevel,
      uint_t                                     maxLevel,
      uint_t                                     maxKrylowDim    = 1,
      uint_t                                     restartLength   = 1,
      real_t                                     arnoldiTOL      = 1e-16,
      real_t                                     approxTOL       = 1e-16,
      real_t                                     doubleOrthoTOL  = 0,
      std::shared_ptr< Solver< OperatorType > >  preconditioner  = std::make_shared< IdentityPreconditioner< OperatorType > >()
      )
      : storage_( storage )
      , minLevel_( minLevel )
      , maxLevel_( maxLevel )
      , maxKrylowDim_( maxKrylowDim )
      , restartLength_( restartLength )
      , arnoldiTOL_( arnoldiTOL )
      , approxTOL_( approxTOL )
      , doubleOrthoTOL_( doubleOrthoTOL )
      , preconditioner_( preconditioner )
      , flag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
      , printInfo_( false )
      , r0_( "r0", storage_, minLevel_, maxLevel_ )
      , rPrec_( "rPrec", storage_, minLevel_, maxLevel_ )
      , wPrec_( "wPrec", storage_, minLevel_, maxLevel_ )
      , orthoDiff_( "orthoDiff", storage_, minLevel_, maxLevel_ )
      , timingTree_( storage->getTimingTree() )
    {}

    ~GMRESSolver() = default;

    // GMRES ACCORDING TO Y.SAAD

    void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
    {
        timingTree_->start( "GMRES Solver" );

        double approxERR = approxTOL_ + 1;
        bool isOnFirstRestart = true;

        init(A, x, b, level, true);

        for(uint_t j = 1; j < maxKrylowDim_; j++)
        {
            if(j % restartLength_ == 0 )
            {
                generateFinalApproximation( x, level );
                init(A, x, b, level, false );
                isOnFirstRestart = false;
                WALBERLA_LOG_INFO_ON_ROOT( " [GMRES] restarted " );
                continue;
            }

            int currentIndex = j % restartLength_;

            if(isOnFirstRestart)
            {
                FunctionType w( "w", storage_, minLevel_, maxLevel_ );
                w.copyBoundaryConditionFromFunction( x );
                vecV_.push_back( w );
            } else
            {
                vecV_[currentIndex].interpolate( 0.0 , level, flag_ );
            }
            A.apply( vecV_[currentIndex - 1], vecV_[currentIndex], level, flag_ );
            wPrec_.interpolate( 0.0, level, hyteg::Inner );
            wPrec_.copyBoundaryConditionFromFunction( x );
            preconditioner_->solve( A, wPrec_, vecV_[currentIndex], level );
            vecV_[currentIndex].assign( {1.0}, {wPrec_}, level, flag_ );
            H_ = matrixResize( H_, currentIndex + 1, currentIndex );
            for(uint_t i = 1; i <= currentIndex; i++)
            {
                double h = vecV_[currentIndex].dotGlobal( vecV_[i - 1], level, flag_ );
                H_(i-1, currentIndex - 1) = h;
                vecV_[currentIndex].add( {-h}, {vecV_[i-1]}, level, flag_ );
            }
            A.apply( vecV_[currentIndex - 1], orthoDiff_, level, flag_ );
            orthoDiff_.add( {-1.0}, {vecV_[currentIndex]}, level, flag_ );
            if( std::sqrt(orthoDiff_.dotGlobal( orthoDiff_, level, flag_ )) > doubleOrthoTOL_ )
            {
                if( printInfo_ ) {
                    WALBERLA_LOG_INFO_ON_ROOT(" [GMRES] invoked double-orthogonalization at iteration " << j );
                }
                 for(uint_t i = 1; i <= (currentIndex); i++)
                {
                    double h = vecV_[currentIndex].dotGlobal( vecV_[i-1], level, flag_ );
                    H_(i-1, currentIndex - 1) += h;
                    vecV_[currentIndex].add( {-h}, {vecV_[i % restartLength_ - 1]}, level, flag_ );
                }
            }

            double wNorm = std::sqrt( vecV_[currentIndex].dotGlobal( vecV_[currentIndex], level, flag_ ));
            H_(currentIndex, currentIndex - 1) = wNorm;
            y_ = hessenbergMinimizer( beta_, H_, Q_, 1, approxERR );
            if( printInfo_ ) {
                WALBERLA_LOG_INFO_ON_ROOT(" [GMRES] approximated residual after " << j << " iterations : " << approxERR );
            }
            vecV_[currentIndex].assign( {1.0/wNorm}, {vecV_[currentIndex]}, level, flag_ );

            if(wNorm <= arnoldiTOL_ || approxERR <= approxTOL_)
            {
                break;
            }
        }
        generateFinalApproximation( x, level );

        timingTree_->stop( "GMRES Solver" );
        return;
    }

private:
    void init(const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level, bool isFirstInit)
    {
        A.apply(x, r0_, level, flag_);
        r0_.assign( {1.0, -1.0}, {b, r0_}, level, flag_ );
        rPrec_.interpolate( 0.0, level, flag_ );
        preconditioner_->solve( A, rPrec_, r0_, level );
        r0_.assign( {1.0}, {rPrec_}, level, flag_ );

        r0_.copyBoundaryConditionFromFunction( x );
        rPrec_.copyBoundaryConditionFromFunction( x );
        wPrec_.copyBoundaryConditionFromFunction( x );
        orthoDiff_.copyBoundaryConditionFromFunction( x );

        beta_ = std::sqrt( r0_.dotGlobal( r0_, level, flag_ ) );

        H_ = Eigen::MatrixXd::Zero(0,0);
        Q_ = Eigen::MatrixXd::Ones(1,1);

        if(isFirstInit)
        {
            FunctionType v0( "v0", storage_, minLevel_, maxLevel_ );
            v0.assign( {1.0/beta_}, {r0_}, level, flag_ );
            vecV_.push_back( v0 );
        } else
        {
            vecV_[0].assign( {1.0/beta_}, {r0_}, level, flag_ );
        }
    }

    const std::shared_ptr< PrimitiveStorage >& storage_;
    uint_t minLevel_;
    uint_t maxLevel_;
    uint_t maxKrylowDim_;
    uint_t restartLength_;
    real_t arnoldiTOL_;
    real_t approxTOL_;
    real_t doubleOrthoTOL_;
    hyteg::DoFType flag_;
    bool printInfo_;
    std::shared_ptr< Solver< OperatorType > > preconditioner_;

    double beta_;
    std::vector< FunctionType > vecV_;
    Eigen::MatrixXd H_;
    Eigen::MatrixXd Q_;
    Eigen::VectorXd y_;
    FunctionType r0_;
    FunctionType rPrec_;
    FunctionType wPrec_;
    FunctionType orthoDiff_;

    std::shared_ptr< walberla::WcTimingTree > timingTree_;

    void generateFinalApproximation(const FunctionType& x, uint_t level)
    {
        for(int i = 0; i < y_.size(); i++)
        {
            x.add( {y_(i)}, {vecV_[i]}, level, flag_ );
        }
    }

    Eigen::VectorXd getUnitVector(int length, int j)
    {
        Eigen::VectorXd answer = Eigen::VectorXd::Zero(length);
        answer(j) = 1;
        return answer;
    }

    Eigen::VectorXd getUnitVector(int length, int j, double val)
    {
        Eigen::VectorXd answer = Eigen::VectorXd::Zero(length);
        answer(j) = val;
        return answer;
    }

   Eigen::MatrixXd getIdentity(int rowsNew, int colsNew) {
        Eigen::MatrixXd answer = Eigen::MatrixXd::Zero(rowsNew, colsNew);
        for(int i = 0; i < std::min(rowsNew, colsNew); i++)
        {
            answer(i, i) = 1;
        }
        return answer;
    }

    Eigen::MatrixXd matrixResize(Eigen::MatrixXd inputMatrix, int rowsNew, int colsNew)
    {
      Eigen::MatrixXd outputMatrix = Eigen::MatrixXd::Zero(rowsNew, colsNew);
        for(uint_t j = 0; j < std::min(rowsNew, (int) inputMatrix.rows()); j++)
        {
            for(int i = 0; i < std::min(colsNew, (int) inputMatrix.cols()); i++)
            {
                outputMatrix(j, i) = inputMatrix(j, i);
            }
        }
      return outputMatrix;
    }

    Eigen::MatrixXd expandQMatrix(Eigen::MatrixXd inputMatrix, int expandBy)
    {
        Eigen::MatrixXd outputMatrix = getIdentity(inputMatrix.rows() + expandBy, inputMatrix.cols() + expandBy);
        for(int i = 0; i < inputMatrix.rows(); i++)
        {
            for(int j = 0; j < inputMatrix.cols(); j++)
            {
                outputMatrix(i,j) = inputMatrix(i,j);
            }
        }
        return outputMatrix;
    }

    Eigen::VectorXd triangSolver(Eigen::MatrixXd triMat, Eigen::VectorXd targetVector)
    {
        Eigen::VectorXd answer = Eigen::VectorXd::Zero(triMat.cols());
        for(int i = triMat.cols()-1; i >= 0; i--)
        {
            double targetOffset = 0;
            for(int j = triMat.cols()-1; j > i; j--)
            {
                targetOffset += triMat(i, j) * answer(j);
            }
            answer(i) = (targetVector(i) - targetOffset) / triMat(i,i);
        }
        return answer;
    }

    Eigen::VectorXd hessenbergMinimizer(double beta, Eigen::MatrixXd& H, Eigen::MatrixXd& Q, int numUnfinishedColumns, double& approxERR)
    {
        Eigen::VectorXd approxVector;

        if(H.rows() != Q.rows())
        {
            Q = expandQMatrix(Q, H.rows() - Q.rows());
        }

        for(int i = H.cols() - numUnfinishedColumns; i < H.cols(); i++)
        {
            H.col(i) = Q * H.col(i);
        }

        for(int i = H.cols() - numUnfinishedColumns; i < H.cols(); i++)
        {
            if(H(i+1, i) == 0)
            {
                continue;
            }
            double s = H(i+1, i) / std::sqrt( std::pow(H(i, i), 2) + std::pow(H(i+1, i), 2) );
            double c = H(i, i) / std::sqrt( std::pow(H(i, i), 2) + std::pow(H(i+1, i), 2) );
            Eigen::MatrixXd QNew = getIdentity( H.rows(), H.rows() );
            QNew(i, i) = c;
            QNew(i+1, i+1) = c;
            QNew(i+1, i) = -s;
            QNew(i, i+1) = s;

            H = QNew * H;
            Q = QNew * Q;
        }

        Eigen::MatrixXd equationLeftSide = H.block(0, 0, H.cols(), H.cols());
        Eigen::VectorXd targetVector = Q * getUnitVector(H.rows(), 0, beta);
        Eigen::VectorXd equationRightSide = targetVector.head(H.cols());
        approxERR = std::abs(targetVector(targetVector.rows() - 1));

        return triangSolver(equationLeftSide, equationRightSide);
    }
};

} // namespace hyteg
