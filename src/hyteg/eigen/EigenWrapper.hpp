/*
 * Copyright (c) 2017-2019 Daniel Drzisga.
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

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/HytegDefinitions.hpp"


namespace walberla {
namespace mpi {

/// Serialisation of Eigen::DenseMatrix for MPI communication
template < typename T, // Element type of SendBuffer
           typename G, // Growth policy of SendBuffer
           typename EigenScalarType,
           int numRows,
           int numCols,
           int options >
GenericSendBuffer< T, G >& operator<<( GenericSendBuffer< T, G >&                                buffer,
                                       const Eigen::Matrix<EigenScalarType, numRows, numCols, options> &eigenMatrix)
{
   for ( int rowIdx = 0; rowIdx < numRows; ++rowIdx )
   {
      for ( int colIdx = 0; colIdx < numCols; ++colIdx )
      {
         buffer << eigenMatrix( rowIdx, colIdx );
      }
   }
   return buffer;
}

/// Deserialisation of Eigen::DenseMatrix for MPI communication
template < typename T, // Element type of RecvBuffer
           typename EigenScalarType,
           int numRows,
           int numCols,
           int options >
GenericRecvBuffer< T >& operator>>( GenericRecvBuffer< T >&                             buffer,
                                    Eigen::Matrix<EigenScalarType, numRows, numCols, options> eigenMatrix) {
   for ( int rowIdx = 0; rowIdx < numRows; ++rowIdx )
   {
      for ( int colIdx = 0; colIdx < numCols; ++colIdx )
      {
         buffer >> eigenMatrix( rowIdx, colIdx );
      }
   }
   return buffer;
}

} // namespace mpi
} // namespace walberla
