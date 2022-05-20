/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include "hyteg/indexing/CouplingCount.hpp"

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/edgedofspace/EdgeDoFPetsc.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/indexing/CouplingCountFreeFunction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

// This test checks whether we obtain the same number of non-zero matrix entries,
// when we assemble the operator as a matrix in PETSc or when we ask the operator
// object to count its own couplings.

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using namespace hyteg;

template < typename opType, template < class > class fType >
void compareCounts( std::shared_ptr< PrimitiveStorage > storage, std::string tag, const uint_t maxLevel, bool beVerbose )
{
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( "Testing " << tag );
   WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------" );

   PETScManager petscManager;

   for ( uint_t level = 0; level <= maxLevel; level++ )
   {
      // determine indices and dimensions
      fType< idx_t > enumerator( "enumerator", storage, level, level );
      enumerator.enumerate( level );
      typedef typename FunctionTrait< fType< idx_t > >::Tag enumTag;

      uint_t globalDoFs = numberOfGlobalDoFs< enumTag >( *storage, level );

      if ( beVerbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "LEVEL = " << level );
         WALBERLA_LOG_INFO_ON_ROOT( "Number of DoFs = " << globalDoFs );
      }

      // setup operator
      opType oper( storage, level, level );
      uint_t nnzHyTeG = indexing::getNumberOfGlobalDoFCouplings( oper, level );

      // assemble matrix
      PETScSparseMatrix< opType > petscMat;
      petscMat.createMatrixFromOperator( oper, level, enumerator, All );
      uint_t nnzPETSc = petscMat.getInfo().getNNZ();

      if ( beVerbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "HyTeG: #couplings/nnz ....... " << nnzHyTeG );
         WALBERLA_LOG_INFO_ON_ROOT( "PETSc: #couplings/nnz ....... " << nnzPETSc );
      }
      WALBERLA_ASSERT_EQUAL( nnzHyTeG, nnzPETSc );
   }
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // -------------------
   //  General behaviour
   // -------------------
   bool   beVerbose  = true;
   bool   run2DTests = true;
   bool   run3DTests = true;
   uint_t maxLevel   = 3;

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   if ( run2DTests )
   {
      // std::string           meshFileName = "../../data/meshes/quad_16el.msh";
      // std::string           meshFileName = "../../data/meshes/tri_1el.msh";
      std::string           meshFileName = "../../data/meshes/annulus_coarse.msh";
      MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
      SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      loadbalancing::roundRobin( setupStorage );
      std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      if ( beVerbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------" );
         WALBERLA_LOG_INFO_ON_ROOT( " Primitives for 2D Tests                                " );
         WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------" );
         WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );
      }

      // run the tests
      compareCounts< P1ConstantMassOperator, P1Function >( storage, "P1-P1 Scalar Operator (2D mesh)", maxLevel, beVerbose );
      compareCounts< P2ConstantMassOperator, P2Function >( storage, "P2-P2 Scalar Operator (2D mesh)", maxLevel, beVerbose );
      compareCounts< P2P1TaylorHoodStokesOperator, P2P1TaylorHoodFunction >(
          storage, "P2-P1 Taylor Hood Stokes (2D mesh)", maxLevel, beVerbose );

      typedef EdgeDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
          P2EdgeDoFMassOperator;
      compareCounts< P2EdgeDoFMassOperator, EdgeDoFFunction >( storage, "EdgeDoF-EdgeDoF (2D mesh)", maxLevel, beVerbose );
   }

   // ----------------------------
   //  Prepare setup for 3D tests
   // ----------------------------
   if ( run3DTests )
   {
      // std::string           meshFileName = "../../data/meshes/3D/tet_tilted_1el.msh";
      // std::string           meshFileName = "../../data/meshes/3D/pyramid_2el.msh";
      // std::string           meshFileName = "../../data/meshes/3D/three_tets_with_two_joint_faces.msh";
      // std::string           meshFileName = "../../data/meshes/3D/pyramid_4el.msh";
      std::string meshFileName = "../../data/meshes/3D/regular_octahedron_8el.msh";
      // std::string           meshFileName = "../../data/meshes/3D/cube_4120el.msh";
      MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
      SetupPrimitiveStorage setupStorage3D( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      setupStorage3D.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
      loadbalancing::roundRobin( setupStorage3D );
      std::shared_ptr< PrimitiveStorage > storage3D = std::make_shared< PrimitiveStorage >( setupStorage3D );

      if ( beVerbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------" );
         WALBERLA_LOG_INFO_ON_ROOT( " Primitives for 3D Tests                                " );
         WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------" );
         WALBERLA_LOG_INFO_ON_ROOT( "Meshfile generated from '" << meshFileName << "'" );
         WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage3D );
      }

      // We start by simply using a P1-P1 operator
      compareCounts< P1ConstantMassOperator, P1Function >( storage3D, "P1-P1 (3D mesh)", maxLevel, beVerbose );
      typedef EdgeDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
          P2EdgeDoFMassOperator;
      compareCounts< P2EdgeDoFMassOperator, EdgeDoFFunction >( storage3D, "E-E (3D mesh)", maxLevel, beVerbose );
      compareCounts< P2ConstantMassOperator, P2Function >( storage3D, "P2-P2 (3D mesh)", maxLevel, beVerbose );
      compareCounts< P2P1TaylorHoodStokesOperator, P2P1TaylorHoodFunction >(
          storage3D, "P2-P1 Taylor Hood Stokes (3D mesh)", maxLevel, beVerbose );
   }

   return 0;
}
