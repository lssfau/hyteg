/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "EGFunction.hpp"

#include "hyteg/sparseassembly/DirichletBCs.hpp"

namespace hyteg {

void applyDirichletBC( const EGFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   //WALBERLA_LOG_INFO( "Warning EDG applies dirichlet boundary conditions for P1!"
   //                   "This is could (should?) be done weakly with Nietsche's method." )

   applyDirichletBC( numerator.getConformingPart()->component( 0 ), mat, level );
   applyDirichletBC( numerator.getConformingPart()->component( 1 ), mat, level );
   if ( numerator.getStorage()->hasGlobalCells() )
      applyDirichletBC( numerator.getConformingPart()->component( 2 ), mat, level );
}

} // namespace hyteg
