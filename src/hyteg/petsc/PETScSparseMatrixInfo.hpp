/*
 * Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
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

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixInfo.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

class PETScSparseMatrixInfo : public SparseMatrixInfo
{
 public:
   PETScSparseMatrixInfo() = delete;

   PETScSparseMatrixInfo( MatInfo info )
   : info_( info )
   {
      PETScManager::ensureIsInitialized();
   };

   /// returns the number of (used) non-zero entries in the matrix
   uint_t getNNZ() const final { return info_.nz_used; };

   friend std::ostream& operator<<( std::ostream& os, const PETScSparseMatrixInfo& matInfo );

 private:
   MatInfo info_;
};

inline std::ostream& operator<<( std::ostream& os, const PETScSparseMatrixInfo& matInfo )
{
   os << "* block size ............................. " << matInfo.info_.block_size << "\n"
      << "* number of nonzeros (alloced) ........... " << matInfo.info_.nz_allocated << "\n"
      << "* number of nonzeros (used) .............. " << matInfo.info_.nz_used << "\n"
      << "* number of nonzeros (unneeded) .......... " << matInfo.info_.nz_unneeded << "\n"
      << "* memory allocated ....................... " << matInfo.info_.memory << "\n"
      << "* no. of matrix assemblies called ........ " << matInfo.info_.assemblies << "\n"
      << "* no. of mallocs during MatSetValues() ... " << matInfo.info_.mallocs << "\n";
   return os;
}

} // namespace hyteg

#endif
