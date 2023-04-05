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

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

#include "array"
namespace hyteg {
using walberla::realIsEqual;

bool sortVertices( Point3D& iCoord, Point3D& nCoord )
{
   if ( walberla::realIsEqual( iCoord[0], nCoord[0] ) )
   {
      if ( walberla::realIsEqual( iCoord[1], nCoord[1] ) )
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
   if ( realIsEqual(first[0][0],second[0][0]) && realIsEqual(first[0][1],second[0][1]) && realIsEqual(first[0][2],second[0][2]))
   {
      if ( realIsEqual(first[1][0],second[1][0]) && realIsEqual(first[1][1],second[1][1]) && realIsEqual(first[1][2],second[1][2]))
      {
         if ( realIsEqual(first[2][0],second[2][0]) && realIsEqual(first[2][1],second[2][1]) && realIsEqual(first[2][2],second[2][2]))
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

void runtest( uint_t coarseRefinements, uint_t level, const std::string& meshFile )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Running: coarseRef:" << coarseRefinements << "Level:" << level << "Mesh:" << meshFile )
   auto originalMeshInfo =
       meshFile == "TORUS" ?
           hyteg::MeshInfo::meshTorus(
               4, 3, real_c( 6.2 ), { 3 }, real_c(0.0), real_c( 2.0 * walberla::math::pi / walberla::real_c( 2 * 6 ) ) ) :
           hyteg::MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage6el( originalMeshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage6el = std::make_shared< PrimitiveStorage >( setupStorage6el );

   auto                  meshInfoRefined = MeshInfo::refinedCoarseMesh( originalMeshInfo, coarseRefinements );
   SetupPrimitiveStorage setupStorageCoarsen( meshInfoRefined, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storageRefined = std::make_shared< PrimitiveStorage >( setupStorageCoarsen );

   //celldof::allCellTypes
   auto normalVertexVector  = createVertexVector( storage6el, level );
   auto refinedVertexVector = createVertexVector( storageRefined, level - coarseRefinements );
   WALBERLA_CHECK_EQUAL( normalVertexVector.size(), refinedVertexVector.size() )

   for ( uint_t i = 0; i < normalVertexVector.size(); i++ )
   {
      for ( uint_t n = 0; n < 4; n++ )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( toVec3( normalVertexVector[i][n] ),
                               toVec3( refinedVertexVector[i][n] ),
                               normalVertexVector[i][0] << " " << normalVertexVector[i][1] << " " << normalVertexVector[i][2]
                                                        << " " << normalVertexVector[i][3] << "\n"
                                                        << refinedVertexVector[i][0] << " " << refinedVertexVector[i][1] << " "
                                                        << refinedVertexVector[i][2] << " " << refinedVertexVector[i][3] )
      }
   }
}
} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::runtest( 1, 1, "../../data/meshes/3D/cube_6el.msh" );
   hyteg::runtest( 1, 2, "../../data/meshes/3D/cube_6el.msh" );
   hyteg::runtest( 1, 3, "../../data/meshes/3D/cube_6el.msh" );
   hyteg::runtest( 2, 2, "../../data/meshes/3D/cube_6el.msh" );
   hyteg::runtest( 1, 2, "TORUS" );
}
