/*
 * Copyright (c) 2020-2023 Marcus Mohr, Daniel Bauer.
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

#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.hpp"

#include "core/DataTypes.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using namespace hyteg;

void printTestHdr( std::string msg )
{
   std::string sep = "------------------------------------------------------";
   WALBERLA_LOG_INFO_ON_ROOT( "" << sep << "\n" << msg << ":" );
}

template < class op_T, class rowSumForm_T, bool isP1 >
struct factory
{};

template < class op_T, class rowSumForm_T >
struct factory< op_T, rowSumForm_T, true >
{
   static op_T genOperator( std::shared_ptr< PrimitiveStorage >& storage, uint_t level, rowSumForm_T& form )
   {
      WALBERLA_UNUSED( form );
      return op_T( storage, level, level );
   };
};

template < class op_T, class rowSumForm_T >
struct factory< op_T, rowSumForm_T, false >
{
   static op_T genOperator( std::shared_ptr< PrimitiveStorage >& storage, uint_t level, rowSumForm_T& form )
   {
      return op_T( storage, level, level, form );
   };
};

template < class cOperType, class vOperType, class RowSumFormType, bool isP1 >
void compareOperators( std::shared_ptr< PrimitiveStorage >& storage,
                       uint_t                               level,
                       std::shared_ptr< RowSumFormType >&   rowSumForm,
                       real_t                               bound,
                       bool                                 outputVTK = false )
{
   cOperType cOper = factory< cOperType, RowSumFormType, isP1 >::genOperator( storage, level, *rowSumForm );
   vOperType vOper( storage, level, level, rowSumForm );

   typedef typename cOperType::srcType funcType;

   funcType funcInp( "input", storage, level, level );
   funcType funcOut1( "output 1", storage, level, level );
   funcType funcOut2( "output 2", storage, level, level );
   funcType funcErr( "error", storage, level, level );

   funcInp.interpolate( 1.0, level, All );

   cOper.apply( funcInp, funcOut1, level, All );
   vOper.apply( funcInp, funcOut2, level, All );

   funcErr.assign( { 1.0, -1.0 }, { funcOut1, funcOut2 }, level, All );
   real_t maxErr = funcErr.getMaxMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "--> Maximal difference = " << maxErr );

   if ( outputVTK )
   {
      VTKOutput vtkOutput( "../../output", "DiagonalOperatorTest", storage );
      vtkOutput.add( funcInp );
      vtkOutput.add( funcOut1 );
      vtkOutput.add( funcOut2 );
      vtkOutput.add( funcErr );
      vtkOutput.write( level );
   }

   WALBERLA_CHECK_LESS( maxErr, bound );
}

#ifdef HYTEG_BUILD_WITH_PETSC

template < class operType, class RowSumFormType >
void testAssembly( std::shared_ptr< PrimitiveStorage >& storage, uint_t level, std::shared_ptr< RowSumFormType >& rowSumForm )
{
   PETScManager                  petscManager;
   PETScSparseMatrix< operType > matrix( "diagonal matrix" );

   typename operType::srcType::template FunctionType< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   operType oper( storage, level, level, rowSumForm );
   matrix.createMatrixFromOperator( oper, level, enumerator, All );

   WALBERLA_LOG_INFO_ON_ROOT( "--> Sparse matrix assembly worked" );
}

template < class cOperType, class vOperType, class RowSumFormType, bool isP1 >
void compareMatrices( std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                               level,
                      std::shared_ptr< RowSumFormType >&   rowSumForm,
                      real_t                               bound )
{
   PETScManager                   petscManager;
   PETScSparseMatrix< vOperType > testMat( "diagonal matrix 1" );
   PETScSparseMatrix< cOperType > compMat( "diagonal matrix 2" );

   typename vOperType::srcType::template FunctionType< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   vOperType vOper( storage, level, level, rowSumForm );
   testMat.createMatrixFromOperator( vOper, level, enumerator, All );

   cOperType cOper = factory< cOperType, RowSumFormType, isP1 >::genOperator( storage, level, *rowSumForm );
   compMat.createMatrixFromOperator( cOper, level, enumerator, All );

   uint_t nDiagVals = numberOfGlobalDoFs< typename vOperType::srcType::Tag >( *storage, level );

   MatInfo infoTest;
   MatGetInfo( testMat.get(), MAT_GLOBAL_SUM, &infoTest );
   uint_t nDiagTest = uint_c( infoTest.nz_used );
   WALBERLA_CHECK_EQUAL( nDiagTest, nDiagVals );

   // cannot compare the two, as the constant operator does not produce a true diagonal matrix
   // MatInfo infoComp;
   // MatGetInfo( compMat.get(), MAT_GLOBAL_SUM, &infoComp );
   // uint_t nDiagComp = uint_c( infoComp.nz_used );
   // WALBERLA_ASSERT_EQUAL( nDiagTest, nDiagComp );

   bool beVerbose = true;
   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      MatInfo info;
      MatGetInfo( testMat.get(), MAT_GLOBAL_SUM, &info );
      WALBERLA_LOG_INFO_ON_ROOT( "Info on testMat:" );
      WALBERLA_LOG_INFO_ON_ROOT( "* block size ............................. " << real_c( info.block_size ) );
      WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros ..................... " << info.nz_used );
      WALBERLA_LOG_INFO_ON_ROOT( "* memory allocated ....................... " << info.memory );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of matrix assemblies called ........ " << info.assemblies );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of mallocs during MatSetValues() ... " << info.mallocs << "\n" );

      MatGetInfo( compMat.get(), MAT_GLOBAL_SUM, &info );
      WALBERLA_LOG_INFO_ON_ROOT( "Info on compMat:" );
      WALBERLA_LOG_INFO_ON_ROOT( "* block size ............................. " << real_c( info.block_size ) );
      WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros ..................... " << info.nz_used << " (no true diagonal matrix)" );
      WALBERLA_LOG_INFO_ON_ROOT( "* memory allocated ....................... " << info.memory );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of matrix assemblies called ........ " << info.assemblies );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of mallocs during MatSetValues() ... " << info.mallocs << "\n" );
   }

   // determine difference between matrices and its norms
   PetscErrorCode ierr;
   ierr = MatAXPY( compMat.get(), -1.0, testMat.get(), DIFFERENT_NONZERO_PATTERN );
   if ( ierr != 0 )
   {
      WALBERLA_ABORT( "Shit happened in PETSc! Our fault most likely!" );
   }

   PetscReal normFrb = 0.0;
   MatNorm( compMat.get(), NORM_FROBENIUS, &normFrb );

   PetscReal normOne = 0.0;
   MatNorm( compMat.get(), NORM_1, &normOne );

   PetscReal normInf = 0.0;
   MatNorm( compMat.get(), NORM_INFINITY, &normInf );

   if ( beVerbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Norms of difference matrix:" );
      WALBERLA_LOG_INFO_ON_ROOT( "* Frobenius norm ...... " << normFrb );
      WALBERLA_LOG_INFO_ON_ROOT( "* 1-norm .............. " << normOne );
      WALBERLA_LOG_INFO_ON_ROOT( "* Infinity norm ....... " << normInf << "\n" );
   }

   // export operators for diagnosis
   // exportOperator< cOperType >( cOper, "DiagonalMatrix_cOper.m", "cMat", storage, level, false, false, true );
   // exportOperator< vOperType >( vOper, "DiagonalMatrix_vOper.m", "vMat", storage, level, false, false, true );

   std::array< real_t, 3 > limits = { bound, bound, bound };

   WALBERLA_CHECK_LESS_EQUAL( normFrb, limits[0] );
   WALBERLA_CHECK_LESS_EQUAL( normOne, limits[1] );
   WALBERLA_CHECK_LESS_EQUAL( normInf, limits[2] );
}

template < class cOperType, class vOperType, class FormType, bool isP1 >
void compareDiagonals( std::shared_ptr< PrimitiveStorage >& storage,
                       uint_t                               level,
                       std::shared_ptr< FormType >&         form,
                       real_t                               bound )
{
   PETScManager                   petscManager;
   PETScSparseMatrix< vOperType > testMat( "diagonal matrix 1" );
   PETScSparseMatrix< cOperType > compMat( "diagonal matrix 2" );

   typename vOperType::srcType::template FunctionType< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   vOperType vOper( storage, level, level, form );
   testMat.createMatrixFromOperator( vOper, level, enumerator, All );

   cOperType cOper = factory< cOperType, FormType, isP1 >::genOperator( storage, level, *form );
   compMat.createMatrixFromOperator( cOper, level, enumerator, All );

   PetscInt localSize, globalSize;
   MatGetSize( testMat.get(), &localSize, &globalSize );
   MatSetSizes( testMat.get(), localSize, localSize, globalSize, globalSize );
   MatSetSizes( compMat.get(), localSize, localSize, globalSize, globalSize );

   Vec testDiag;
   VecCreate( walberla::mpi::MPIManager::instance()->comm(), &testDiag );
   VecSetType( testDiag, VECMPI );
   VecSetSizes( testDiag, localSize, globalSize );
   VecSetUp( testDiag );
   MatGetDiagonal( testMat.get(), testDiag );

   Vec compDiag;
   VecCreate( walberla::mpi::MPIManager::instance()->comm(), &compDiag );
   VecSetType( compDiag, VECMPI );
   VecSetSizes( compDiag, localSize, globalSize );
   VecSetUp( compDiag );
   MatGetDiagonal( compMat.get(), compDiag );

   // determine difference between diagonals and its norms
   PetscErrorCode ierr;
   ierr = VecAXPY( compDiag, -1.0, testDiag );
   if ( ierr != 0 )
   {
      WALBERLA_ABORT( "Shit happened in PETSc! Our fault most likely!" );
   }

   PetscReal normOne = 0.0;
   VecNorm( compDiag, NORM_1, &normOne );

   PetscReal normTwo = 0.0;
   VecNorm( compDiag, NORM_2, &normTwo );

   PetscReal normInf = 0.0;
   VecNorm( compDiag, NORM_INFINITY, &normInf );

   WALBERLA_LOG_INFO_ON_ROOT( "Norms of difference vector:" );
   WALBERLA_LOG_INFO_ON_ROOT( "* 1-norm .............. " << normOne );
   WALBERLA_LOG_INFO_ON_ROOT( "* 2-norm .............. " << normTwo );
   WALBERLA_LOG_INFO_ON_ROOT( "* Infinity norm ....... " << normInf << "\n" );

   WALBERLA_CHECK_LESS_EQUAL( normInf, bound );
}

#endif

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // --------------------------
   //  Prepare underlying forms
   // --------------------------
   auto p1MassFormFenics =
       std::make_shared< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >();
   std::shared_ptr< P1RowSumForm > lumpedMassFormP1 = std::make_shared< P1RowSumForm >( p1MassFormFenics );

   auto p2MassFormFenics =
       std::make_shared< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >();
   std::shared_ptr< P2RowSumForm > lumpedMassFormP2 = std::make_shared< P2RowSumForm >( p2MassFormFenics );

   auto                            p1MassFormHyTeG3D       = std::make_shared< forms::p1_mass_blending_q4 >();
   std::shared_ptr< P1RowSumForm > lumpedMassFormP1HyTeG3D = std::make_shared< P1RowSumForm >( p1MassFormHyTeG3D );

   auto                            p2MassFormHyTeG       = std::make_shared< forms::p2_mass_blending_q4 >();
   std::shared_ptr< P2RowSumForm > lumpedMassFormP2HyTeG = std::make_shared< P2RowSumForm >( p2MassFormHyTeG );

   typedef P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > P1LaplaceForm_T;
   auto p1LaplaceForm2D = std::make_shared< P1LaplaceForm_T >();

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   std::string           meshFileName = "../../data/meshes/quad_16el.msh";
   MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t level = 2;

   // -----------------------
   //  Perform 2D experiment
   // -----------------------

   WALBERLA_LOG_INFO_ON_ROOT( "============\n  2D TESTS\n============" );

   printTestHdr( "Testing Mass Lumping for P1" );
   compareOperators< P1LumpedMassOperator, P1BlendingLumpedDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, real_c( std::is_same< real_t, double >() ? 5e-17 : 1e-9 ) );

   printTestHdr( "Testing Inverted Mass Lumping for P1" );
   compareOperators< P1LumpedInvMassOperator, P1BlendingLumpedInverseDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, real_c( std::is_same< real_t, double >() ? 1e-10 : 5e-5 ), true );

   printTestHdr( "Testing Mass Lumping for P2" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingLumpedDiagonalOperator, P2RowSumForm, false >(
       storage, level, lumpedMassFormP2, real_c( 1e-17 ) );

   printTestHdr( "Testing Laplace for P1" );
   typedef DiagonalNonConstantOperator< P1ElementwiseOperator, P1LaplaceForm_T, false > P1BlendingLaplaceDiagonalOperator;
   compareOperators< P1DiagonalLaplaceOperator, P1BlendingLaplaceDiagonalOperator, P1LaplaceForm_T, true >(
       storage, level, p1LaplaceForm2D, real_c( 5e-15 ) );

   // -------------------------------------
   //  Perform 2D experiment with blending
   // -------------------------------------
   // Regression test for https://i10git.cs.fau.de/hyteg/hyteg/-/merge_requests/665

#ifdef HYTEG_BUILD_WITH_PETSC

   WALBERLA_LOG_INFO_ON_ROOT( "=======================\n  2D TESTS (blending)\n=======================" );

   meshInfo = MeshInfo::meshAnnulus( 2, 4, MeshInfo::CRISS, 6, 2 );
   SetupPrimitiveStorage setupStorageBlending( meshInfo,
                                               walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   AnnulusMap::setMap( setupStorageBlending );
   std::shared_ptr< PrimitiveStorage > storageBlending = std::make_shared< PrimitiveStorage >( setupStorageBlending );

   printTestHdr( "Testing Mass Diagonal for P2" );
   std::shared_ptr< forms::p2_mass_blending_q5 > p2MassFormBlending = std::make_shared< forms::p2_mass_blending_q5 >();
   compareDiagonals< P2ElementwiseOperator< forms::p2_mass_blending_q5 >,
                     DiagonalNonConstantOperator< P2ElementwiseOperator, forms::p2_mass_blending_q5 >,
                     forms::p2_mass_blending_q5,
                     false >(
       storageBlending, level, p2MassFormBlending, real_c( std::is_same< real_t, double >() ? 8e-18 : 5e-09 ) );

#endif

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

   WALBERLA_LOG_INFO_ON_ROOT( "============\n  3D TESTS\n============" );

   printTestHdr( "Testing Mass Lumping for P1 (FEniCS Form)" );
   compareOperators< P1LumpedMassOperator, P1BlendingLumpedDiagonalOperator, P1RowSumForm, true >(
       storage3D, level, lumpedMassFormP1, real_c( std::is_same< real_t, double >() ? 1e-16 : 5e-09 ) );

   printTestHdr( "Testing Inverted Mass Lumping for P1 (FEniCS Form)" );
   compareOperators< P1LumpedInvMassOperator, P1BlendingLumpedInverseDiagonalOperator, P1RowSumForm, true >(
       storage3D, level, lumpedMassFormP1, real_c( std::is_same< real_t, double >() ? 2e-11 : 5e-4 ) );

   printTestHdr( "Testing Mass Lumping for P2 (FEniCS Form)" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingLumpedDiagonalOperator, P2RowSumForm, false >(
       storage3D, level, lumpedMassFormP2, real_c( std::is_same< real_t, double >() ? 2e-18 : 4e-7 ) );

   printTestHdr( "Testing Mass Lumping for P1 (HyTeG Form)" );
   compareOperators< P1LumpedMassOperator, P1BlendingLumpedDiagonalOperator, P1RowSumForm, true >(
       storage3D, level, lumpedMassFormP1HyTeG3D, real_c( std::is_same< real_t, double >() ? 1e-16 : 6e-9 ) );

   printTestHdr( "Testing Mass Lumping for P2 (HyTeG Form)" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingLumpedDiagonalOperator, P2RowSumForm, false >(
       storage3D, level, lumpedMassFormP2HyTeG, real_c( std::is_same< real_t, double >() ? 4e-17 : 6e-9 ) );

   // ----------------------
   //  Test Matrix Assembly
   // ----------------------

#ifdef HYTEG_BUILD_WITH_PETSC

   WALBERLA_LOG_INFO_ON_ROOT( "===================\n  Matrix Assembly\n===================" );

   printTestHdr( "Testing Mass Lumping for P2 (HyTeG Form, 2D)" );
   testAssembly< P2BlendingLumpedDiagonalOperator, P2RowSumForm >( storage, level, lumpedMassFormP2HyTeG );

   printTestHdr( "Testing Mass Lumping for P1 (FEniCS Form, 2D)" );
   compareMatrices< P1LumpedMassOperator, P1BlendingLumpedDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, real_c( std::is_same< real_t, double >() ? 1e-16 : 2e-08 ) );

   printTestHdr( "Testing Inverted Mass Lumping for P1 (FEniCS Form, 2D)" );
   compareMatrices< P1LumpedInvMassOperator, P1BlendingLumpedInverseDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, real_c( std::is_same< real_t, double >() ? 1e-11 : 2e-04 ) );

#endif

   return 0;
}
