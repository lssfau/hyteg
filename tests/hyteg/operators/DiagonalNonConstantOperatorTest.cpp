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

#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.hpp"

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/FunctionTraits.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/ElementwiseOperatorPetsc.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
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
   // cOperType cOper( storage, level, level );
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

   auto p1PSPGFormFenics =
       std::make_shared< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >();
   std::shared_ptr< P1RowSumForm > pspgFormP1 = std::make_shared< P1RowSumForm >( p1PSPGFormFenics );

   auto p1MassFormHyTeG3D = std::make_shared< P1Form_mass3D >();
   std::shared_ptr< P1RowSumForm > lumpedMassFormP1HyTeG3D = std::make_shared< P1RowSumForm >( p1MassFormHyTeG3D );

   auto p2MassFormHyTeG = std::make_shared< P2Form_mass >();
   std::shared_ptr< P2RowSumForm > lumpedMassFormP2HyTeG = std::make_shared< P2RowSumForm >( p2MassFormHyTeG );

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   std::string           meshFileName = "../../data/meshes/quad_16el.msh";
   MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t level = 4;

   // -----------------------
   //  Perform 2D experiment
   // -----------------------

   WALBERLA_LOG_INFO_ON_ROOT( "============\n  2D TESTS\n============" );

   printTestHdr( "Testing Mass Lumping for P1" );
   compareOperators< P1LumpedMassOperator, P1BlendingDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-17 );

   printTestHdr( "Testing Inverted Mass Lumping for P1" );
   compareOperators< P1LumpedInvMassOperator, P1BlendingInverseDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-10, true );

   // doesn't make sense numerically, but just to check another operator different from mass
   printTestHdr( "Testing PSPG for P1" );
   compareOperators< P1LumpedMassOperator, P1BlendingDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-17 );

   printTestHdr( "Testing Mass Lumping for P2" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingDiagonalOperator, P2RowSumForm, false >(
       storage, level, lumpedMassFormP2, 1e-17 );

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

   printTestHdr( "Testing Mass Lumping for P1" );
   compareOperators< P1LumpedMassOperator, P1BlendingDiagonalOperator, P1RowSumForm, true >(
       storage3D, level, lumpedMassFormP1, 1e-17 );

   printTestHdr( "Testing Inverted Mass Lumping for P1" );
   compareOperators< P1LumpedInvMassOperator, P1BlendingInverseDiagonalOperator, P1RowSumForm, true >(
       storage, level, lumpedMassFormP1, 1e-10 );

   printTestHdr( "Testing Mass Lumping for P2" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingDiagonalOperator, P2RowSumForm, false >(
       storage, level, lumpedMassFormP2, 1e-17 );

   printTestHdr( "Testing Mass Lumping for P1 (HyTeG Form)" );
   compareOperators< P1LumpedMassOperator, P1BlendingDiagonalOperator, P1RowSumForm, true >(
       storage3D, level, lumpedMassFormP1HyTeG3D, 1e-17 );

   printTestHdr( "Testing Mass Lumping for P2 (HyTeG Form)" );
   compareOperators< P2ConstantRowSumOperator, P2BlendingDiagonalOperator, P2RowSumForm, false >(
       storage, level, lumpedMassFormP2HyTeG, 1e-17 );

   return 0;
}
