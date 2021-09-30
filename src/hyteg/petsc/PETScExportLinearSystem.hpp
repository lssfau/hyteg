/*
 * Copyright (c) 2017-2019 Nils Kohl, Dominik Thoennes, Marcus Mohr.
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

#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

/// \brief Exports the operator and right-hand side of a linear system.
///
/// Uses PETSc functionality to generate the matrix associated with the given operator
/// and the right-hand side vector and export both to a file in ASCII format which can
/// be executed as a Matlab function.
///
/// \param op                 the operator to export
/// \param rhs                the right-hand side function
/// \param dirichletSolution  a function that has the solution interpolated on the Dirichlet boundary
/// \param fileMatrix         name of the file to write the matrix to
/// \param nameMatrix         name of the matrix, must be different from base part of fileMatrix;
///                           when the script is executed in Matlab this will be the final
///                           name of the sparse matrix
/// \param fileRHS            name of the file to write the right-hand side vector to
/// \param nameRHS            name of the right-hand side, must be different from base part of fileRHS;
/// \param storage            primitive storage
/// \param level              refinement level to be considered
/// \param elimDirichletBC    whether to zero row/columns for Dirichlet boundary values,
///                           if true, this also adapts the right-hand side vector
/// \param beVerbose          should function be talkative or not
///
template < class OperatorType, template < class > class FunctionType, class FunctionTag >
void exportLinearSystem( OperatorType                        op,
                         const FunctionType< real_t >&       rhs,
                         const FunctionType< real_t >&       dirichletSolution,
                         std::string                         fileMatrix,
                         std::string                         nameMatrix,
                         std::string                         fileRHS,
                         std::string                         nameRHS,
                         std::shared_ptr< PrimitiveStorage > storage,
                         uint_t                              level,
                         bool                                elimDirichletBC,
                         bool                                beVerbose = false,
                         bool                                binary    = false,
                         PetscViewerFormat                   format    = PETSC_VIEWER_ASCII_MATRIXMARKET )
{
   // Get dimension of function space
   uint_t localDoFs  = hyteg::numberOfLocalDoFs< FunctionTag >( *storage, level );
   uint_t globalDoFs = hyteg::numberOfGlobalDoFs< FunctionTag >( *storage, level );
   if ( localDoFs != globalDoFs )
   {
      WALBERLA_ABORT( "localDoFs and globalDoFs must agree for this app!" );
   }
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Dimension of function space is " << globalDoFs );
   }

   // Fire up PETSc
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Firing up PETSc" );
   }
   PETScManager pmgr;

   // Create PETSc matrix
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Converting Operator to PETSc matrix" );
   }
   PETScSparseMatrix< OperatorType > petscMatrix( localDoFs, globalDoFs, nameMatrix.c_str() );
   FunctionType< idx_t >             numerator( "numerator", storage, level, level );
   numerator.enumerate( level );
   petscMatrix.createMatrixFromOperator( op, level, numerator );

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Converting RHS to PETSc vector" );
   }
   PETScVector< real_t, FunctionType > petscRHS( rhs, numerator, level, All, nameRHS );

   // Zero rows and columns of "Dirichlet DoFs"
   if ( elimDirichletBC )
   {
      petscMatrix.applyDirichletBCSymmetrically( dirichletSolution, numerator, petscRHS, level );
   }

   // Write out matrix
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Exporting Operator to file '" << fileMatrix << "'" );
   }
   petscMatrix.print( fileMatrix.c_str(), binary, format );

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * Exporting RHS to file '" << fileRHS << "'" );
   }
   petscRHS.print( fileRHS.c_str(), binary, format );
}

} // namespace hyteg

#endif
