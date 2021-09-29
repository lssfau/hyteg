/*
 * Copyright (c) 2017-2019 Marcus Mohr.
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


/// \file
/// This source code file contains free functions for eigenvalue estimation.

#pragma once

#include "hyteg/eigen/EigenWrapper.hpp"
#include "core/Abort.h"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/flags.hpp"

namespace hyteg {

  using walberla::real_t;
  using walberla::uint_t;

  // =================================================================================================

  /// \brief Compute an estimate for the spectral radius of an operator using the Power Method.
  ///
  /// The estimate for the spetral radius is computed using a simple Power iteration.
  /// Note that for this to work the initial vector itrVec must not be "deficient", i.e. it
  /// must contain as component a multiple of the eigenvector associated with the eigenvalue
  /// of largest modulus.
  /// For details see e.g. the book "Matrix Computation" by Gene Golub and Charles van Loan.
  ///
  /// \param op          operator to be investigated
  /// \param itrVec      auxilliary vector needed for performing iterations
  /// \param numIts      number of iterations to be performed
  /// \param storage     primitive storage object
  /// \param level       grid level to work on (operator & gridfunctions)
  ///
  /// \return Estimate for the spectral radius of op  
  template < typename OperatorType >
  real_t estimateSpectralRadiusWithPowerIteration( OperatorType& op,
                                                   typename OperatorType::srcType &itrVec,
                                                   typename OperatorType::srcType &auxVec,
                                                   const uint_t numIts,
                                                   const std::shared_ptr<PrimitiveStorage> &storage,
                                                   const uint_t level ) {

    // scale input vector
    real_t norm = std::sqrt( itrVec.dotGlobal( itrVec, level, hyteg::All ) );
    itrVec.assign( {1.0/norm}, {itrVec}, level, hyteg::All );

    // prepare iteration
    op.apply( itrVec, auxVec, level, hyteg::All );
    real_t radius = 0.0;

    // run power iteration
    for( uint_t it = 1; it <= numIts; ++it ) {
      norm = std::sqrt( auxVec.dotGlobal( auxVec, level, hyteg::All ) );
      itrVec.assign( {1.0/norm}, {auxVec}, level, hyteg::All );
      op.apply( itrVec, auxVec, level, hyteg::All );
      radius = itrVec.dotGlobal( auxVec, level, hyteg::All );
    }

    return radius;
  }

  // =================================================================================================

  /// \brief Computes estimates for the smallest and largest eigenvalue of a symmetric
  /// positive definite operator
  ///
  /// The extremal eigenvalues are estimated by setting up the tridiagonal matrix associated
  /// with the Lanczos process. This is done via the conjugate gradient method. The eigenvalues
  /// of this matrix are computed with the help of the Eigen library. While the matrix setup
  /// requires parallel communication, the eigenvalue computation is done by each MPI process
  /// individually and does not involve commmunication.
  /// For algorithmic details see e.g. the book "Iterative methods for sparse linear systems"
  /// by Yousef Saad.
  ///
  /// \param op          operator to be investigated
  /// \param itrVec      auxilliary vector needed for performing CG iterations
  /// \param rhsVec      right-hand side vector used for CG iterations
  /// \param numIts      number of CG steps performed corresponds to dimension of matrix
  /// \param storage     primitive storage object
  /// \param level       grid level to work on (operator & gridfunctions)
  /// \param lowerBound  on return contains an estimate for the smallest eigenvalue of op
  /// \param lowerBound  on return contains an estimate for the largest eigenvalue of op
  template < typename OperatorType >
  void estimateSpectralBoundsWithCG( OperatorType op,
                                     typename OperatorType::srcType &itrVec,
                                     typename OperatorType::srcType &rhsVec,
                                     const uint_t numIts,
                                     const std::shared_ptr<PrimitiveStorage> &storage,
                                     const uint_t level,
                                     real_t& lowerBound,
                                     real_t& upperBound ) {

    // setup CG solver and generate the Lanczos matrix
    CGSolver<OperatorType> cg = CGSolver<OperatorType>( storage, level, level );
    std::vector<real_t> mainDiag, subDiag;
    cg.setupLanczosTriDiagMatrix( op, itrVec, rhsVec, level, numIts, mainDiag, subDiag );

    // compute spectrum of matrix using Eigen library
    Eigen::Map<Eigen::VectorXd> dVec( mainDiag.data(), numIts );
    Eigen::Map<Eigen::VectorXd> sVec( subDiag.data(), numIts - 1 );
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.computeFromTridiagonal( dVec, sVec, Eigen::EigenvaluesOnly );

    // Eigen sorts eigenvalues ascendingly, extract smallest and largest one
    Eigen::VectorXd ev = es.eigenvalues();
    lowerBound = ev[0];
    upperBound = ev[numIts-1];
  }

} // namespace hyteg

