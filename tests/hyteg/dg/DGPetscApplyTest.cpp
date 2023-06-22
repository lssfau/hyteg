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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;

namespace hyteg {

void dgPetscApplyTest( uint_t level, const MeshInfo& meshInfo, real_t eps )
{
   using namespace dg;

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto basis    = std::make_shared< DGBasisLinearLagrange_Example >();
   auto diffForm = std::make_shared< DGDiffusionForm_Example >( ( storage->hasGlobalCells() ? 0.5 : 1 ) );
   auto massForm = std::make_shared< DGMassForm_Example >();

   DGFunction< real_t > src( "src", storage, level, level, basis, 1 );
   DGFunction< real_t > tmp( "tmp", storage, level, level, basis, 1 );
   DGFunction< real_t > hytegDst( "hytegDst", storage, level, level, basis, 1 );
   DGFunction< real_t > petscDst( "petscDst", storage, level, level, basis, 1 );
   DGFunction< real_t > err( "error", storage, level, level, basis, 1 );
   DGFunction< idx_t >  numerator( "numerator", storage, level, level, basis, 1 );

   DGOperator L( storage, level, level, diffForm );
   DGOperator M( storage, level, level, massForm );

   std::function< real_t( const hyteg::Point3D& ) > srcFunction = []( const hyteg::Point3D& x ) {
      return x[0] * x[0] * x[0] * x[0] * std::sinh( x[1] ) * std::cos( x[2] );
   };

   // interpolation of some function into the src vector
   tmp.evaluateLinearFunctional( srcFunction, level );
   PETScCGSolver< DGOperator > solverM( storage, level, numerator );
   solverM.solve( M, src, tmp, level );

   numerator.enumerate( level );

   const uint_t globalDoFs = src.getNumberOfGlobalDoFs( level );
   WALBERLA_LOG_INFO_ON_ROOT( "Global DoFs: " << globalDoFs );

   // HyTeG apply
   L.apply( src, hytegDst, level, All, Replace );

   // PETSc apply
   PETScVector< real_t, DGFunction > srcPetscVec;
   PETScVector< real_t, DGFunction > dstPetscVec;
   PETScSparseMatrix< DGOperator >   petscMatrix;

   srcPetscVec.createVectorFromFunction( src, numerator, level );
   dstPetscVec.createVectorFromFunction( petscDst, numerator, level );
   petscMatrix.createMatrixFromOperator( L, level, numerator );

   WALBERLA_CHECK( petscMatrix.isSymmetric() );

   MatMult( petscMatrix.get(), srcPetscVec.get(), dstPetscVec.get() );

   dstPetscVec.createFunctionFromVector( petscDst, numerator, level );

   // compare
   err.assign( { 1.0, -1.0 }, { hytegDst, petscDst }, level );
   auto maxMag = err.getMaxMagnitude( level );

   WALBERLA_LOG_INFO_ON_ROOT( "Error max mag = " << maxMag );

   // VTK
   //   VTKOutput vtkOutput( "../../output", "DGPetscApplyTest", storage );
   //   vtkOutput.add( src );
   //   vtkOutput.add( hytegDst );
   //   vtkOutput.add( petscDst );
   //   vtkOutput.add( err );
   //   vtkOutput.write( level, 0 );

   WALBERLA_CHECK_LESS( maxMag, eps );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   hyteg::dgPetscApplyTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 6.0e-15 );
   hyteg::dgPetscApplyTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/annulus_coarse.msh" ), 8.0e-14 );
   hyteg::dgPetscApplyTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-16 );
   hyteg::dgPetscApplyTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ), 1.0e-15 );
   hyteg::dgPetscApplyTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ), 1.0e-15 );
   hyteg::dgPetscApplyTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/regular_octahedron_8el.msh" ), 1.0e-15 );

   return EXIT_SUCCESS;
}
