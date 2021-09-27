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
#include "hyteg/elementwiseoperators/ElementwiseOperatorPetsc.hpp"

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

// This test sets up the matrix associated with a constant operator
// in sparse format using the interface to PETSc and compares it to
// the matrix we get for the corresponding elementwise operator.
// We do that for several operators/forms in 2D and 3D.

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using namespace hyteg;

template < typename OpTypeConst, typename OpTypeElem, template < class > class FuncType >
void compareMatrices( std::shared_ptr< PrimitiveStorage > storage,
                      std::string                         tag,
                      const uint_t                        level,
                      std::array< real_t, 3 >             limits,
                      bool                                beVerbose )
{
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------" );
   WALBERLA_LOG_INFO_ON_ROOT( " Running matrix comparison for: " << tag );
   WALBERLA_LOG_INFO_ON_ROOT( "-------------------------------------------------" );

   // determine indices and dimensions
   FuncType< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   typedef typename FunctionTrait< FuncType< idx_t > >::Tag    enumTag;
   uint_t                                                      globalDoFs = numberOfGlobalDoFs< enumTag >( *storage, level );
   uint_t                                                      localDoFs  = numberOfLocalDoFs< enumTag >( *storage, level );

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Matrices will have dimension " << globalDoFs << "\n" );
   }

   // setup operators
   OpTypeElem  elemWiseOp( storage, level, level );
   OpTypeConst constantOp( storage, level, level );

   // assemble matrices
   PETScSparseMatrix< OpTypeElem, FuncType > elemWisePETScMat( localDoFs, globalDoFs );
   elemWisePETScMat.createMatrixFromOperator( elemWiseOp, level, enumerator, All );

   PETScSparseMatrix< OpTypeConst, FuncType > constantPETScMat( localDoFs, globalDoFs );
   constantPETScMat.createMatrixFromOperator( constantOp, level, enumerator, All );

   // determine difference between matrices and their norms
   // PetscErrorCode MatAXPY(Mat Y,PetscReal a,Mat X,MatStructure str)  // SAME_NONZERO_PATTERN
   PetscErrorCode ierr;
   ierr = MatAXPY( constantPETScMat.get(), -1.0, elemWisePETScMat.get(), DIFFERENT_NONZERO_PATTERN );
   if ( ierr != 0 )
   {
      WALBERLA_ABORT( "Shit happened in PETSc! Our fault most likely!" );
   }

   PetscReal normFrb = 0.0;
   MatNorm( constantPETScMat.get(), NORM_FROBENIUS, &normFrb );

   PetscReal normOne = 0.0;
   MatNorm( constantPETScMat.get(), NORM_1, &normOne );

   PetscReal normInf = 0.0;
   MatNorm( constantPETScMat.get(), NORM_INFINITY, &normInf );

   if ( beVerbose )
   {
      // MatInfo info;
      WALBERLA_LOG_INFO_ON_ROOT( "Info on constantPETScMat:" );
      WALBERLA_LOG_INFO_ON_ROOT( "" << constantPETScMat.getInfo() );
      // MatGetInfo( constantPETScMat.get(), MAT_GLOBAL_SUM, &info );
      // WALBERLA_LOG_INFO_ON_ROOT( "Info on constantPETScMat:" );
      // WALBERLA_LOG_INFO_ON_ROOT( "* block size ............................. " << real_c( info.block_size ) );
      // WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros (alloced) ........... " << info.nz_allocated );
      // WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros (used) .............. " << info.nz_used      );
      // WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros (unneeded) .......... " << info.nz_unneeded  );
      // WALBERLA_LOG_INFO_ON_ROOT( "* memory allocated ....................... " << info.memory );
      // WALBERLA_LOG_INFO_ON_ROOT( "* no. of matrix assemblies called ........ " << info.assemblies );
      // WALBERLA_LOG_INFO_ON_ROOT( "* no. of mallocs during MatSetValues() ... " << info.mallocs << "\n" );

      WALBERLA_LOG_INFO_ON_ROOT( "Info on elemWisePETScMat:" );
      WALBERLA_LOG_INFO_ON_ROOT( "" << elemWisePETScMat.getInfo() );
      // MatGetInfo( elemWisePETScMat.get(), MAT_GLOBAL_SUM, &info );
      // WALBERLA_LOG_INFO_ON_ROOT( "Info on elemWisePETScMat:" );
      // WALBERLA_LOG_INFO_ON_ROOT( "* block size ............................. " << real_c( info.block_size ) );
      // WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros (alloced) ........... " << info.nz_allocated );
      // WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros (used) .............. " << info.nz_used      );
      // WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros (unneeded) .......... " << info.nz_unneeded  );
      // WALBERLA_LOG_INFO_ON_ROOT( "* memory allocated ....................... " << info.memory );
      // WALBERLA_LOG_INFO_ON_ROOT( "* no. of matrix assemblies called ........ " << info.assemblies );
      // WALBERLA_LOG_INFO_ON_ROOT( "* no. of mallocs during MatSetValues() ... " << info.mallocs << "\n" );

      WALBERLA_LOG_INFO_ON_ROOT( "Norms of difference matrix:" );
      WALBERLA_LOG_INFO_ON_ROOT( "* Frobenius norm ...... " << normFrb );
      WALBERLA_LOG_INFO_ON_ROOT( "* 1-norm .............. " << normOne );
      WALBERLA_LOG_INFO_ON_ROOT( "* Infinity norm ....... " << normInf << "\n" );
   }

   WALBERLA_CHECK_LESS_EQUAL( normFrb, limits[0] );
   WALBERLA_CHECK_LESS_EQUAL( normOne, limits[1] );
   WALBERLA_CHECK_LESS_EQUAL( normInf, limits[2] );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );
   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   std::string           meshFileName = "../../data/meshes/quad_16el.msh";
   MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   bool beVerbose = false;
   // bool beVerbose = true;
   uint_t level = 4;

   // -------------------
   //  Run some 2D tests
   // -------------------
   compareMatrices< P1ConstantLaplaceOperator, P1ElementwiseLaplaceOperator, P1Function >(
       storage, "P1Laplace", level, {1e-12, 1e-13, 1e-13}, beVerbose );
   compareMatrices< P1ConstantMassOperator, P1ElementwiseMassOperator, P1Function >(
       storage, "P1Mass", level, {1e-16, 1e-17, 1e-17}, beVerbose );

   level = 3;
   compareMatrices< P2ConstantLaplaceOperator, P2ElementwiseLaplaceOperator, P2Function >(
       storage, "P2Laplace", level, {1.5e-12, 1.5e-13, 1.5e-13}, beVerbose );
   compareMatrices< P2ConstantMassOperator, P2ElementwiseMassOperator, P2Function >(
       storage, "P2Mass", level, {1e-16, 1e-17, 1e-17}, beVerbose );
   compareMatrices< P2ConstantLaplaceOperator, P2ElementwiseDivKGradOperator, P2Function >(
       storage, "P2DivKGrad", level, {1.5e-12, 1e-13, 1.5e-13}, beVerbose );
   compareMatrices< P2P1TaylorHoodStokesOperator, P2P1ElementwiseConstantCoefficientStokesOperator, P2P1TaylorHoodFunction >(
       storage, "P2P1StokesCC", level, {2.0e-12, 1.5e-13, 1.5e-13}, beVerbose );

   // ----------------------------
   //  Prepare setup for 3D tests
   // ----------------------------
   meshFileName                     = "../../data/meshes/3D/pyramid_tilted_4el.msh";
   MeshInfo              meshInfo3D = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage3D( meshInfo3D, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage3D );
   std::shared_ptr< PrimitiveStorage > storage3D = std::make_shared< PrimitiveStorage >( setupStorage3D );

   // -------------------
   //  Run some 3D tests
   // -------------------
   level = 3;
   compareMatrices< P1ConstantLaplaceOperator, P1ElementwiseLaplaceOperator, P1Function >(
       storage3D, "P1Laplace - 3D", level, {1e-12, 1e-13, 1e-13}, beVerbose );
   compareMatrices< P1ConstantMassOperator, P1ElementwiseMassOperator, P1Function >(
       storage3D, "P1Mass - 3D", level, {1e-16, 1e-17, 1e-17}, beVerbose );

   compareMatrices< P2ConstantLaplaceOperator, P2ElementwiseLaplaceOperator, P2Function >(
       storage3D, "P2Laplace - 3D", level, {1e-13, 1e-13, 1e-13}, beVerbose );
   compareMatrices< P2ConstantMassOperator, P2ElementwiseMassOperator, P2Function >(
       storage3D, "P2Mass - 3D", level, {1e-17, 1e-18, 1e-18}, beVerbose );
   compareMatrices< P2P1TaylorHoodStokesOperator, P2P1ElementwiseConstantCoefficientStokesOperator, P2P1TaylorHoodFunction >(
       storage3D, "P2P1StokesCC", level, {1.5e-12, 1.5e-13, 1.5e-13}, beVerbose );

   return 0;
}
