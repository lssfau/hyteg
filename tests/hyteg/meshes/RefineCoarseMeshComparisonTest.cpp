/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Nils Kohl.
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

#include "core/debug/CheckFunctions.h"
#include "core/math/Constants.h"

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "array"
namespace hyteg {
bool sortVertices( Point3D& iCoord, Point3D& nCoord )
{
   if ( walberla::floatIsEqual( iCoord[0], nCoord[0] ) )
   {
      if ( walberla::floatIsEqual( iCoord[1], nCoord[1] ) )
      {
         return iCoord[2] < nCoord[2];
      }
      else
      {
         return iCoord[1] < nCoord[1];
      }
   }
   else
   {
      return iCoord[0] < nCoord[0];
   }
}

bool sortVerticesArray( std::array< Point3D, 4 >& first, std::array< Point3D, 4 >& second )
{
   if ( toVec3( first[0] ) == toVec3( second[0] ) )
   {
      if ( toVec3( first[1] ) == toVec3( second[1] ) )
      {
         if ( toVec3( first[2] ) == toVec3( second[2] ) )
         {
            if ( sortVertices( first[3], second[3] ) )
            {
               return true;
            }
            else
            {
               return false;
            }
         }
         else
         {
            if ( sortVertices( first[2], second[2] ) )
            {
               return true;
            }
            else
            {
               return false;
            }
         }
      }
      else
      {
         if ( sortVertices( first[1], second[1] ) )
         {
            return true;
         }
         else
         {
            return false;
         }
      }
   }
   else
   {
      if ( sortVertices( first[0], second[0] ) )
      {
         return true;
      }
      else
      {
         return false;
      }
   }
}

std::vector< std::array< Point3D, 4 > > createVertexVector( std::shared_ptr< PrimitiveStorage >& storage, uint_t level )
{
   std::vector< std::array< Point3D, 4 > > normalVertexVector;
   for ( const auto& cell : storage->getCells() )
   {
      for ( const auto& cType : celldof::allCellTypes )
      {
         for ( const auto& micro : celldof::macrocell::Iterator( level, cType, 0 ) )
         {
            std::array< Point3D, 4 > newVertices;
            uint_t                   idx               = 0;
            auto                     microCellVertices = celldof::macrocell::getMicroVerticesFromMicroCell( micro, cType );
            for ( indexing::Index& microVertexIndex : microCellVertices )
            {
               newVertices[idx] = vertexdof::macrocell::coordinateFromIndex( level, *cell.second, microVertexIndex );
               idx++;
            }
            std::sort( newVertices.begin(), newVertices.end(), sortVertices );
            normalVertexVector.emplace_back( newVertices );
         }
      }
   }
   std::sort( normalVertexVector.begin(), normalVertexVector.end(), sortVerticesArray );
   return normalVertexVector;
}

void runtest( uint_t coarseRefinements, uint_t level, const MeshInfo& originalMeshInfo )
{
   SetupPrimitiveStorage setupStorage6el( originalMeshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage6el = std::make_shared< PrimitiveStorage >( setupStorage6el );

   auto                  meshInfoCoarsen = MeshInfo::refinedCoarseMesh( originalMeshInfo, coarseRefinements );
   SetupPrimitiveStorage setupStorageCoarsen( meshInfoCoarsen, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storageCoarsen = std::make_shared< PrimitiveStorage >( setupStorageCoarsen );

   //celldof::allCellTypes
   auto normalVertexVector  = createVertexVector( storage6el, level );
   auto coarsenVertexVector = createVertexVector( storageCoarsen, level - coarseRefinements );
   WALBERLA_CHECK_EQUAL( normalVertexVector.size(), coarsenVertexVector.size() )

   for ( uint_t i = 0; i < normalVertexVector.size(); i++ )
   {
      for ( uint_t n = 0; n < 4; n++ )
      {
         WALBERLA_CHECK_EQUAL( toVec3( normalVertexVector[i][n] ),
                               toVec3( coarsenVertexVector[i][n] ),
                               normalVertexVector[i][0] << " " << normalVertexVector[i][1] << " " << normalVertexVector[i][2]
                                                        << " " << normalVertexVector[i][3] << "\n"
                                                        << coarsenVertexVector[i][0] << " " << coarsenVertexVector[i][1] << " "
                                                        << coarsenVertexVector[i][2] << " " << coarsenVertexVector[i][3] )
      }
   }
}
} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::runtest( 1, 1, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ) );
   hyteg::runtest( 1, 2, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ) );
   hyteg::runtest( 1, 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ) );
   hyteg::runtest( 2, 2, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" ) );
   hyteg::runtest( 1,
                   2,
                   hyteg::MeshInfo::meshTorus(
                       4, 3, 6.2, { 3 }, 0.0, 2.0 * walberla::math::pi / walberla::real_c( 2 * 6 ) ) );
}