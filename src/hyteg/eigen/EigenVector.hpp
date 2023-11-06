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

#include "core/mpi/MPIManager.h"

#include "hyteg/eigen/EigenVectorProxy.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

/// \brief Builds and returns an Eigen VectorX from a HyTeG FE function.
template < template < class > class FunctionType >
VectorXr createEigenVectorFromFunction( const FunctionType< real_t >& src,
                                        const FunctionType< idx_t >&  numerator,
                                        uint_t                        level,
                                        DoFType                       flag = All )
{
   WALBERLA_CHECK_EQUAL(
       walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Eigen vectors are not suited for MPI parallel applications." )

   const uint_t rows = numberOfLocalDoFs( numerator, level );

   VectorXr x;
   x.resize( rows );

   auto proxy = std::make_shared< EigenVectorProxy >( x );
   src.toVector( numerator, proxy, level, flag );

   return x;
}

/// \brief Builds and returns an Eigen VectorX from a HyTeG FE function.
template < template < class > class FunctionType >
void assignEigenVectorFromFunction( const FunctionType< real_t >& src,
                                    VectorXr&                     dst,
                                    const FunctionType< idx_t >&  numerator,
                                    uint_t                        level,
                                    DoFType                       flag = All )
{
   WALBERLA_CHECK_EQUAL(
       walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Eigen vectors are not suited for MPI parallel applications." )

   const uint_t rows = numberOfLocalDoFs( numerator, level );
   WALBERLA_CHECK_EQUAL( dst.rows(), rows );

   auto proxy = std::make_shared< EigenVectorProxy >( dst );
   src.toVector( numerator, proxy, level, flag );
}

/// \brief Writes data from an Eigen VectorX back to the coefficients of a HyTeG FE function.
template < template < class > class FunctionType >
void assignFunctionFromEigenVector( const VectorXr&               src,
                                    const FunctionType< real_t >& dst,
                                    const FunctionType< idx_t >&  numerator,
                                    uint_t                        level,
                                    DoFType                       flag = All )
{
   WALBERLA_CHECK_EQUAL(
       walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Eigen vectors are not suited for MPI parallel applications." )

   const uint_t rows = numberOfLocalDoFs( numerator, level );

   WALBERLA_CHECK_EQUAL( rows, src.size(), "Source vector (Eigen's VectorX) must be of same size as target FE function." );

   auto proxy = std::make_shared< EigenConstVectorProxy >( src );
   dst.fromVector( numerator, proxy, level, flag );
}

} // namespace hyteg
