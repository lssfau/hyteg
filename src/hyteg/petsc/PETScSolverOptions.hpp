/*
 * Copyright (c) 2024 Andreas Burkhart
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

#include <memory>

#include "hyteg/solvers/Solver.hpp"

#include "PETScSparseMatrix.hpp"
#include "PETScVector.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

class PETScSolverOptions
{
 public:
   PETScSolverOptions()
   {
      kspType_ = "fgmres";
      pcType_  = "none";
   }

   PETScSolverOptions( std::vector< std::pair< std::string, std::string > >& options )
   : PETScSolverOptions()
   {
      for ( auto& p : options )
      {
         addOption( p.first, p.second );
      }
   }

   template < class ValueType >
   void addOption( std::string PETScOption, ValueType PETScOptionValue )
   {
      while ( PETScOption[0] == '-' )
      {
         PETScOption.erase( 0, 1 );
      }

      std::stringstream PETScOptionValueStr;
      PETScOptionValueStr << PETScOptionValue;

      options_.push_back( { PETScOption, PETScOptionValueStr.str() } );
   }

   PetscErrorCode applyOptions( KSP& ksp, std::string prefix = "" ) const
   {
      for ( auto& p : options_ )
      {
         std::stringstream PETScOption;
         PETScOption << "-" << prefix << p.first;

         PetscCall( PetscOptionsSetValue( NULL, PETScOption.str().c_str(), p.second.c_str() ) );
      }

      return PETSC_SUCCESS;
   }

   void        setKspType( std::string kspType ) { kspType_ = kspType; }
   void        setPcType( std::string pcType ) { pcType_ = pcType; }
   std::string getKspType() const { return kspType_; }
   std::string getPcType() const { return pcType_; }

 private:
   std::vector< std::pair< std::string, std::string > > options_;
   std::string                                          kspType_;
   std::string                                          pcType_;
};

} // namespace hyteg

#endif
