/*
* Copyright (c) 2017-2025 Marcus Mohr.
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

#include "hyteg/dg1functionspace/DG1Function.hpp"

#include <utility>

#include "core/DataTypes.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/Reduce.h"

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/BlendingHelpers.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::uint_t;

template < typename value_t >
void DG1Function< value_t >::interpolate( const std::function< value_t( const Point3D& ) >& expr,
                                          uint_t                                            level,
                                          DoFType                                           dofType ) const
{
   if ( dofType != All && dofType != Inner )
   {
      WALBERLA_ABORT( "DG1Function::interpolate() requires dofType to be 'All' or 'Inner'" );
   }

   if ( std::dynamic_pointer_cast< DGBasisLinearLagrange_Example >( dgFunction_->basis() ) == nullptr )
   {
      WALBERLA_ABORT( "DG1Function::interpolate(): Ooops, not your standard nodal basis!" );
   }

   using ElementNeighborInfo = volumedofspace::indexing::ElementNeighborInfo;

   if ( this->storage_->hasGlobalCells() )
   {
      for ( auto& it : this->getStorage()->getCells() )
      {
         const auto cellID = it.first;

         const auto memLayout = dgFunction_->volumeDoFFunction()->memoryLayout();
         auto       dofs      = dgFunction_->volumeDoFFunction()->dofMemory( cellID, level );
         uint_t     numDofs   = 4;

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               ElementNeighborInfo neighborInfo =
                   ElementNeighborInfo( idxIt, cellType, level, this->getBoundaryCondition(), cellID, this->getStorage() );
               const std::vector< indexing::Index >& vIndex = neighborInfo.elementVertexIndices();

               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const Point3D point = micromesh::microVertexPosition( this->getStorage(), cellID, level, vIndex[i] );
                  dofs[volumedofspace::indexing::index(
                      idxIt.x(), idxIt.y(), idxIt.z(), cellType, i, numDofs, level, memLayout )] = expr( point );
               }
            }
         }
      }
   }
   else
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         const auto faceID = it.first;

         const auto memLayout = dgFunction_->volumeDoFFunction()->memoryLayout();
         auto       dofs      = dgFunction_->volumeDoFFunction()->dofMemory( faceID, level );
         uint_t     numDofs   = 3;

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               ElementNeighborInfo neighborInfo =
                   ElementNeighborInfo( idxIt, faceType, level, this->getBoundaryCondition(), faceID, this->getStorage() );
               const std::vector< indexing::Index >& vIndex = neighborInfo.elementVertexIndices();

               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const Point3D point = micromesh::microVertexPosition( this->getStorage(), faceID, level, vIndex[i] );
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )] =
                      expr( point );
               }
            }
         }
      }
   }
}

/// explicit instantiations
template class DG1Function< double >;
template class DG1Function< float >;
template class DG1Function< int32_t >;
template class DG1Function< int64_t >;

} // namespace hyteg
