/*
* Copyright (c) 2023 Nils Kohl.
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

#include "core/debug/CheckFunctions.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/eigen/EigenSparseMatrixProxy.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

/// \brief Builds and returns an Eigen::SparseMatrix from a HyTeG operator.
template < typename OperatorType >
Eigen::SparseMatrix< real_t, Eigen::RowMajor >
    createEigenSparseMatrixFromOperator( const OperatorType&                                                   op,
                                         uint_t                                                                level,
                                         const typename OperatorType::srcType::template FunctionType< idx_t >& numeratorSrc,
                                         const typename OperatorType::dstType::template FunctionType< idx_t >& numeratorDst,
                                         DoFType                                                               flag = All )
{
   WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(),
                         1,
                         "Eigen sparse matrices are not suited for MPI parallel applications." );

   const uint_t rows = numberOfLocalDoFs( numeratorDst, level );
   const uint_t cols = numberOfLocalDoFs( numeratorSrc, level );

   auto proxy = std::make_shared< EigenSparseMatrixProxy >( rows, cols );
   op.toMatrix( proxy, numeratorSrc, numeratorDst, level, flag );

   return proxy->getSparseMatrix();
}
} // namespace hyteg
