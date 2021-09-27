/*
 * Copyright (c) 2020 Marcus Mohr.
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
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/ElementwiseOperatorPetsc.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
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

   walberla::math::seedRandomGenerator( 1234 );
   auto rand = []( const Point3D& ) { return walberla::math::realRandom(); };
   funcInp.interpolate( rand, level, All );
   // funcInp.interpolate( 1.9, level, All );

   cOper.apply( funcInp, funcOut1, level, All );
   vOper.apply( funcInp, funcOut2, level, All );

   funcErr.assign( {1.0, -1.0}, {funcOut1, funcOut2}, level, All );
   real_t maxErr = funcErr.getMaxMagnitude( level );
   WALBERLA_LOG_INFO_ON_ROOT( "--> Maximal difference = " << maxErr );
   WALBERLA_CHECK_LESS( maxErr, bound );

   if ( outputVTK )
   {
      VTKOutput vtkOutput( "../../output", "DiagonalOperatorTest", storage );
      vtkOutput.add( funcInp );
      vtkOutput.add( funcOut1 );
      vtkOutput.add( funcOut2 );
      vtkOutput.add( funcErr );
      vtkOutput.write( level );
   }
}

#ifdef HYTEG_BUILD_WITH_PETSC

template < class operType, class RowSumFormType >
void testAssembly( std::shared_ptr< PrimitiveStorage >& storage, uint_t level, std::shared_ptr< RowSumFormType >& rowSumForm )
{
   PETScManager                                                            petscManager;
   PETScSparseMatrix< operType, operType::srcType::template FunctionType > matrix( storage, level, "diagonal matrix" );

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
   PETScManager                                                              petscManager;
   PETScSparseMatrix< vOperType, vOperType::srcType::template FunctionType > testMat( storage, level, "diagonal matrix 1" );
   PETScSparseMatrix< cOperType, vOperType::srcType::template FunctionType > compMat( storage, level, "diagonal matrix 2" );

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

   std::array< real_t, 3 > limits = {bound, bound, bound};

   WALBERLA_CHECK_LESS_EQUAL( normFrb, limits[0] );
   WALBERLA_CHECK_LESS_EQUAL( normOne, limits[1] );
   WALBERLA_CHECK_LESS_EQUAL( normInf, limits[2] );
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

   auto                            p1MassFormHyTeG3D       = std::make_shared< P1Form_mass3D >();
   std::shared_ptr< P1RowSumForm > lumpedMassFormP1HyTeG3D = std::make_shared< P1RowSumForm >( p1MassFormHyTeG3D );

   auto                            p2MassFormHyTeG       = std::make_shared< P2Form_mass >();
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
       storage, level, lumpedMassFormP1, 5e-17 );

   printTestHdr( "Testing Inverted Mass Lumping for P1" );
   compareOperators< P1LumpedInvMassOperator, P1BlendingLumpedInverseDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-10, true );

   printTestHdr( "Testing Mass Lumping for P2" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingLumpedDiagonalOperator, P2RowSumForm, false >(
       storage, level, lumpedMassFormP2, 1e-17 );

   printTestHdr( "Testing Laplace for P1" );
   typedef DiagonalNonConstantOperator< P1ElementwiseOperator, P1LaplaceForm_T, false > P1BlendingLaplaceDiagonalOperator;
   compareOperators< P1DiagonalLaplaceOperator, P1BlendingLaplaceDiagonalOperator, P1LaplaceForm_T, true >(
       storage, level, p1LaplaceForm2D, 5e-15 );

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
       storage3D, level, lumpedMassFormP1, 1e-16 );

   printTestHdr( "Testing Inverted Mass Lumping for P1 (FEniCS Form)" );
   compareOperators< P1LumpedInvMassOperator, P1BlendingLumpedInverseDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-10 );

   printTestHdr( "Testing Mass Lumping for P2 (FEniCS Form)" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingLumpedDiagonalOperator, P2RowSumForm, false >(
       storage, level, lumpedMassFormP2, 1e-17 );

   printTestHdr( "Testing Mass Lumping for P1 (HyTeG Form)" );
   compareOperators< P1LumpedMassOperator, P1BlendingLumpedDiagonalOperator, P1RowSumForm, true >(
       storage3D, level, lumpedMassFormP1HyTeG3D, 5e-17 );

   printTestHdr( "Testing Mass Lumping for P2 (HyTeG Form)" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingLumpedDiagonalOperator, P2RowSumForm, false >(
       storage, level, lumpedMassFormP2HyTeG, 1e-16 );

   // ----------------------
   //  Test Matrix Assembly
   // ----------------------

#ifdef HYTEG_BUILD_WITH_PETSC

   WALBERLA_LOG_INFO_ON_ROOT( "===================\n  Matrix Assembly\n===================" );

   printTestHdr( "Testing Mass Lumping for P2 (HyTeG Form, 2D)" );
   testAssembly< P2BlendingLumpedDiagonalOperator, P2RowSumForm >( storage, level, lumpedMassFormP2HyTeG );

   printTestHdr( "Testing Mass Lumping for P1 (FEniCS Form, 2D)" );
   compareMatrices< P1LumpedMassOperator, P1BlendingLumpedDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-16 );

   printTestHdr( "Testing Inverted Mass Lumping for P1 (FEniCS Form, 2D)" );
   compareMatrices< P1LumpedInvMassOperator, P1BlendingLumpedInverseDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-11 );

#endif

   return 0;
}
