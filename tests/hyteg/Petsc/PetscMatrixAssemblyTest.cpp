/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1BlendingStokesOperator.hpp"
#include "hyteg/composites/P1EpsilonStokesOperator.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"
#include "hyteg/composites/P1P1UzawaDampingFactorEstimationOperator.hpp"
// #include "hyteg/composites/P1PolynomialBlendingStokesOperator.hpp" < --see issue 159
#include "hyteg/composites/P2P1BlendingTaylorHoodStokesOperator.hpp"
#include "hyteg/composites/P2P1SurrogateTaylorHoodStokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodBlockFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/composites/P2P1UzawaDampingFactorEstimationOperator.hpp"
#include "hyteg/composites/P2P2StabilizedStokesOperator.hpp"
#include "hyteg/composites/P2P2UnstableStokesOperator.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/DiagonalNonConstantOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/p1functionspace/P1SurrogateOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2EpsilonOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2SurrogateOperator.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

/// This test checks whether we can assemble a sparse matrix
/// (for/with PETSc) for the given operators. It is mostly a
/// compile time check.

#ifndef HYTEG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON"
#endif

#define NOT_COMPILING( tag ) WALBERLA_LOG_INFO_ON_ROOT( "* NOT COMPILING: " << tag );

using walberla::real_t;
using namespace hyteg;

void printTestHdr( std::string msg )
{
   std::string sep = "------------------------------------------------------";
   WALBERLA_LOG_INFO_ON_ROOT( "" << sep << "\n" << msg << ":\n" << sep );
}

template < class operType >
void testAssembly( std::shared_ptr< PrimitiveStorage >& storage, uint_t level, std::string tag, bool actuallyTest = true )
{
   if ( actuallyTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * " << tag );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * ~ SKIPPING (BUT COMPILES): " << tag );
      return;
   }

   typename operType::srcType::template FunctionType< idx_t > enumeratorSrc( "enumeratorSrc", storage, level, level );
   typename operType::dstType::template FunctionType< idx_t > enumeratorDst( "enumeratorDst", storage, level, level );
   enumeratorSrc.enumerate( level );
   enumeratorDst.enumerate( level );

   PETScManager                  petscManager;
   PETScSparseMatrix< operType > matrix( tag );

   operType oper( storage, level, level );
   matrix.createMatrixFromOperator( oper, level, enumeratorSrc, enumeratorDst, All );
}

template < template < class > class fKind >
void testAssembly( BlockOperator< fKind< real_t >, fKind< real_t > >& oper,
                   std::shared_ptr< PrimitiveStorage >&               storage,
                   uint_t                                             level,
                   std::string                                        tag,
                   bool                                               actuallyTest = true )
{
   if ( actuallyTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * " << tag );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * ~ SKIPPING (BUT COMPILES): " << tag );
      return;
   }

   typedef BlockOperator< fKind< real_t >, fKind< real_t > > operType;
   fKind< idx_t >                                            enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   PETScManager                  petscManager;
   PETScSparseMatrix< operType > matrix( tag );

   matrix.createMatrixFromOperator( oper, level, enumerator, All );
}

// Version for operators the require a form in their ctors
template < class operType, class formType >
void testAssembly( std::shared_ptr< PrimitiveStorage >& storage,
                   uint_t                               level,
                   std::shared_ptr< formType >&         form,
                   std::string                          tag,
                   bool                                 actuallyTest = true )
{
   if ( actuallyTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * " << tag );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * ~ SKIPPING (BUT COMPILES): " << tag );
      return;
   }

   typename operType::srcType::template FunctionType< idx_t > enumeratorSrc( "enumeratorSrc", storage, level, level );
   typename operType::dstType::template FunctionType< idx_t > enumeratorDst( "enumeratorDst", storage, level, level );
   enumeratorSrc.enumerate( level );
   enumeratorDst.enumerate( level );

   PETScManager                  petscManager;
   PETScSparseMatrix< operType > matrix( tag );

   operType oper( storage, level, level, form );
   matrix.createMatrixFromOperator( oper, level, enumeratorSrc, enumeratorDst, All );
}

// Version for (some?) Surrogate Operators (need an interpolation level in ctor)
template < class operType >
void testAssembly( std::shared_ptr< PrimitiveStorage >& storage,
                   uint_t                               level,
                   uint_t                               ipLevel,
                   std::string                          tag,
                   bool                                 actuallyTest = true )
{
   if ( actuallyTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * " << tag );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * ~ SKIPPING (BUT COMPILES): " << tag );
      return;
   }

   PETScManager                  petscManager;
   PETScSparseMatrix< operType > matrix( tag );

   typename operType::srcType::template FunctionType< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   operType oper( storage, level, level, ipLevel );
   matrix.createMatrixFromOperator( oper, level, enumerator, All );
}

// Version for ProjectNormalOperators
template < class operType >
void testAssembly( uint_t level, std::string tag, bool actuallyTest = true )
{
   if ( actuallyTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * " << tag );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * ~ SKIPPING (BUT COMPILES): " << tag );
      return;
   }

   // blended annulus mesh
   auto                  meshInfo = MeshInfo::meshAnnulus( 0.5, 1.0, MeshInfo::CRISS, 6, 6 );
   SetupPrimitiveStorage annulusStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   annulusStorage.setMeshBoundaryFlagsOnBoundary( 3, 0, true );
   AnnulusMap::setMap( annulusStorage );
   const auto storage = std::make_shared< PrimitiveStorage >( annulusStorage );

   // normal computation
   auto normalInterpolant = []( const Point3D& p ) {
      real_t norm = p.norm();
      real_t sign = ( norm > 0.75 ) ? 1.0 : -1.0;
      return sign / norm * p;
   };

   auto normalFunction = [=]( const Point3D& p, Point3D& n ) -> void { n = normalInterpolant( p ); };

   // standard assemble check
   PETScManager                  petscManager;
   PETScSparseMatrix< operType > matrix( tag );

   typename operType::srcType::template FunctionType< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   operType projectNormalOperator( storage, level, level, normalFunction );
   matrix.createMatrixFromOperator( projectNormalOperator, level, enumerator, All );
}

// -------------------------------
//  Specialised one case versions
// -------------------------------
template <>
void testAssembly< P1ConstantUnsteadyDiffusionOperator >( std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                               level,
                                                          std::string                          tag,
                                                          bool                                 actuallyTest )
{
   if ( actuallyTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * " << tag );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * ~ SKIPPING (BUT COMPILES): " << tag );
      return;
   }

   PETScManager                                             petscManager;
   PETScSparseMatrix< P1ConstantUnsteadyDiffusionOperator > matrix( tag );

   typename P1ConstantUnsteadyDiffusionOperator::srcType::template FunctionType< idx_t > enumerator(
       "enumerator", storage, level, level );
   enumerator.enumerate( level );

   P1ConstantUnsteadyDiffusionOperator oper( storage, level, level, 0.1, 1.0, DiffusionTimeIntegrator::ImplicitEuler );
   matrix.createMatrixFromOperator( oper, level, enumerator, All );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // ----------------------------
   //  Prepare setup for 2D tests
   // ----------------------------
   std::string           meshFileName = "../../data/meshes/quad_16el.msh";
   MeshInfo              meshInfo     = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   uint_t      level = 2;
   std::string rule  = "==============================";

   // ------------------
   //  Scalar Operators
   // ------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n  SCALAR TO SCALAR OPERATORS\n" << rule );

   std::shared_ptr< Form > form = nullptr;

   printTestHdr( "Testing P1 Operators" );

   testAssembly< P1ConstantMassOperator >( storage, level, "P1ConstantOperator" );
   testAssembly< P1SurrogateMassOperator >( storage, level, "P1SurrogateOperator", false );
   testAssembly< P1BlendingMassOperator >( storage, level, "P1VariableOperator", false );

   testAssembly< P1ProjectNormalOperator >( level, "P1ProjectNormalOperator", false );

   printTestHdr( "Testing P2 Operators" );

   testAssembly< P2ConstantMassOperator >( storage, level, "P2ConstantOperator" );

   testAssembly< P2BlendingMassOperator >( storage, level, "P2VariableOperator", false );
   testAssembly< P2SurrogateMassOperator >( storage, level, level, "P2SurrogateOperator", false );
   testAssembly< P2ProjectNormalOperator >( level, "P2ProjectNormalOperator", false );

   // ---------------------
   //  Low-level Operators
   // ---------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n  LOW-LEVEL OPERATORS\n" << rule );
   typedef EdgeDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
       E2EMassOperator;
   testAssembly< E2EMassOperator >( storage, level, "EdgeDoFOperator" );

   typedef EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
       E2VMassOperator;
   testAssembly< E2VMassOperator >( storage, level, "EdgeDoFToVertexDoFOperator" );
   typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
       V2EMassOperator;
   testAssembly< V2EMassOperator >( storage, level, "VertexDoFToEdgeDoFOperator" );

   // -----------------------
   //  Elementwise Operators
   // -----------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n  ELEMENTWISE OPERATORS\n" << rule );

   testAssembly< P1ElementwiseMassOperator >( storage, level, "P1ElementwiseOperator" );
   testAssembly< P2ElementwiseMassOperator >( storage, level, "P2ElementwiseOperator" );

   auto                            p2MassFormHyTeG       = std::make_shared< forms::p2_mass_blending_q4 >();
   std::shared_ptr< P2RowSumForm > lumpedMassFormP2HyTeG = std::make_shared< P2RowSumForm >( p2MassFormHyTeG );
   testAssembly< P2BlendingLumpedDiagonalOperator, P2RowSumForm >(
       storage, level, lumpedMassFormP2HyTeG, "DiagonalNonConstantOperator" );

   testAssembly< P2P1ElementwiseBlendingStokesOperator >( storage, level, "P2P1ElementwiseBlendingStokesOperator" );
   testAssembly< P2ToP1ElementwiseBlendingDivxOperator >( storage, level, "P2ToP1ElementwiseOperator" );
   testAssembly< P1ToP2ElementwiseBlendingDivTxOperator >( storage, level, "P1ToP2ElementwiseOperator" );

   // ---------------------
   //  Composite Operators
   // ---------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n  COMPOSITE OPERATORS\n" << rule );

   // cannot construct operator:
   // testAssembly< P1P1UzawaDampingFactorEstimationOperator >( storage, level, "P1P1UzawaDampingFactorEstimationOperator" );
   // testAssembly< P2P1UzawaDampingFactorEstimationOperator >( storage, level, "P2P1UzawaDampingFactorEstimationOperator" );

   testAssembly< P1ConstantUnsteadyDiffusionOperator >( storage, level, "P1ConstantUnsteadyDiffusionOperator" );

   testAssembly< P1P1StokesOperator >( storage, level, "P1StokesOperator" );
   testAssembly< P2P1TaylorHoodStokesOperator >( storage, level, "P2P1TaylorHoodStokesOperator" );
   testAssembly< P2P2StabilizedStokesOperator >( storage, level, "P2P2StabilizedStokesOperator" );

   testAssembly< P2P1TaylorHoodStokesOperator >( storage, level, "P2P1TaylorHoodStokesOperator" );

   testAssembly< P1BlendingStokesOperator >( storage, level, "P1BlendingStokesOperator", false );
   testAssembly< P1EpsilonStokesOperator >( storage, level, "P1EpsilonStokesOperator" );
   testAssembly< P2P2UnstableStokesOperator >( storage, level, "P2P2UnstableStokesOperator" );
   testAssembly< P2P1BlendingTaylorHoodStokesOperator >( storage, level, "P2P1BlendingTaylorHoodStokesOperator", false );
   testAssembly< P2P1SurrogateTaylorHoodStokesOperator >( storage, level, level, "P2P1SurrogateTaylorHoodStokesOperator", false );

   // ----------------------------
   //  Scalar To Vector Operators
   // ----------------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n  SCALAR TO VECTOR OPERATORS\n" << rule );

   testAssembly< P1ToP2ElementwiseBlendingDivTOperator >( storage, level, "P1ScalarToP2VectorOperator" );
   testAssembly< P1ToP2ConstantDivTxOperator >( storage, level, "P1ToP2ConstantOperator" );

   // ----------------------------
   //  Vector To Scalar Operators
   // ----------------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n  VECTOR TO SCALAR OPERATORS\n" << rule );

   testAssembly< P2ToP1ElementwiseBlendingDivOperator >( storage, level, "P2ScalarToP1VectorOperator" );
   testAssembly< P2ToP1ConstantDivxOperator >( storage, level, "P2ToP1ConstantOperator" );

   // ------------------
   //  Vector Operators
   // ------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n  VECTOR TO VECTOR OPERATORS\n" << rule );
   testAssembly< P1ConstantVectorMassOperator >( storage, level, "P1ConstantVectorMassOperator" );
   testAssembly< P2ConstantVectorMassOperator >( storage, level, "P2ConstantVectorMassOperator" );
   testAssembly< P2ElementwiseBlendingVectorMassOperator >( storage, level, "P2ElementwiseBlendingVectorMassOperator" );

   // WALBERLA_LOG_INFO_ON_ROOT( "Skipping: (since assembly doesn't compile)" );
   // WALBERLA_LOG_INFO_ON_ROOT( " * [all of this kind]" );

   // ------------------
   //  Block Operators
   // ------------------
   WALBERLA_LOG_INFO_ON_ROOT( "" << rule << "\n BLOCK OPERATORS\n" << rule );

   // setup an artificial block operator
   typedef BlockOperator< P2P1TaylorHoodBlockFunction< real_t >, P2P1TaylorHoodBlockFunction< real_t > > myBlockOper;

   myBlockOper oper( storage, level, level, 2, 2 );
   oper.setSubOperator( 0, 0, std::make_shared< OperatorWrapper< P2ConstantVectorLaplaceOperator > >( storage, level, level ) );
   oper.setSubOperator( 1, 0, nullptr );
   oper.setSubOperator( 0, 1, nullptr );
   oper.setSubOperator( 1, 1, std::make_shared< OperatorWrapper< P1ConstantMassOperator > >( storage, level, level ) );

   testAssembly< P2P1TaylorHoodBlockFunction >( oper, storage, level, "artificial diagonal BlockOperator" );

   return 0;
}
