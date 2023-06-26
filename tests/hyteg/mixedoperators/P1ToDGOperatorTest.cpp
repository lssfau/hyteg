/*
* Copyright (c) 2022 Andreas Wagner.
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

#include "hyteg/mixedoperators/P1ToDGOperator.hpp"

#include <hyteg/dg1functionspace/DG1Operator.hpp>
#include <hyteg/dgfunctionspace/DGMassForm_Example.hpp>

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mixedoperators/DGToP1Operator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::math::pi;

void checkP1ToDG1ByIntegral( const uint_t dim )
{
   std::string mesh_file;
   if ( dim == 2 )
      mesh_file = "../../data/meshes/quad_4el.msh";
   else
      mesh_file = "../../data/meshes/3D/cube_6el.msh";

   MeshInfo meshInfo = MeshInfo::fromGmshFile( mesh_file );

   SetupPrimitiveStorage setup( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto            storage = std::make_shared< PrimitiveStorage >( setup, 1 );

   const uint_t level = 2;
   const auto   form  = std::make_shared< P1ToDG1InterpolationForm >();

   P1Function< real_t >  src( "src", storage, level, level );
   DG1Function< real_t > dst( "dst", storage, level, level );
   DG1Function< real_t > Mdst( "Mdst", storage, level, level );

   P1ToDGOperator< P1ToDG1InterpolationForm > op( storage, level, level, form );

   src.interpolate( []( auto p ) { return p[0] - 0.1 * p[1] + 0.2 * p[2]; }, level, All );

   op.apply( src, *dst.getDGFunction(), level, All, hyteg::Replace );

   auto        massForm = std::make_shared< DGMassForm_Example >();
   DG1Operator M( storage, level, level, massForm );

   M.apply( dst, Mdst, level, All, Replace );
   const real_t integralValue         = Mdst.sumGlobal( level, All );
   const real_t expectedIntegralValue = 0.5 - 0.1 * 0.5 + ( dim == 3 ? 0.2 * 0.5 : 0. );

   WALBERLA_CHECK_FLOAT_EQUAL( integralValue, expectedIntegralValue, "integral values must match" );

   // VTK
   // VTKOutput vtkOutput( "../../output", "P1ToDGOperatorTest", storage );
   // vtkOutput.add( src );
   // vtkOutput.add( dst );
   // vtkOutput.write( level, 0 );
}

void enumerateTest()
{
   std::string mesh_file = "../../data/meshes/quad_4el.msh";

   MeshInfo meshInfo = MeshInfo::fromGmshFile( mesh_file );

   SetupPrimitiveStorage setup( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto            storage = std::make_shared< PrimitiveStorage >( setup, 1 );

   const uint_t level = 2;
   const auto   form  = std::make_shared< P1ToDG1InterpolationForm >();

   P1Function< idx_t >  src( "src", storage, level, level );
   DG1Function< idx_t > dst( "dst", storage, level, level );

   P1ToDGOperator< P1ToDG1InterpolationForm, idx_t > op( storage, level, level, form );

   src.enumerate( level );

   op.apply( src, *dst.getDGFunction(), level, All, hyteg::Replace );

   // VTK
   // VTKOutput vtkOutput( "../../output", "P1ToDGOperatorTest", storage );
   // vtkOutput.add( src );
   // vtkOutput.add( dst );
   // vtkOutput.write( level, 0 );
}

void checkTranspose()
{
   std::string mesh_file = "../../data/meshes/quad_4el.msh";

   MeshInfo meshInfo = MeshInfo::fromGmshFile( mesh_file );

   SetupPrimitiveStorage setup( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto            storage = std::make_shared< PrimitiveStorage >( setup, 1 );

   const uint_t level = 2;
   const auto   form  = std::make_shared< P1ToDG1InterpolationForm >();

   P1Function< idx_t >  srcP1( "srcP1", storage, level, level );
   P1Function< idx_t >  dstP1( "dstP1", storage, level, level );
   DG1Function< idx_t > srcDG( "srcDG", storage, level, level );
   DG1Function< idx_t > dstDG( "dstDG", storage, level, level );

   P1ToDGOperator< P1ToDG1InterpolationForm, idx_t > opP1ToDG( storage, level, level, form );
   DGToP1Operator< P1ToDG1InterpolationForm, idx_t > opDGToP1( storage, level, level, form );

   srcP1.enumerate( level );
   srcDG.enumerate( level );

   opP1ToDG.apply( srcP1, *dstDG.getDGFunction(), level, All, hyteg::Replace );
   const auto value1 = dstDG.dotGlobal( srcDG, level, All );

   opDGToP1.apply( *srcDG.getDGFunction(), dstP1, level, All, hyteg::Replace );
   const auto value2 = dstP1.dotGlobal( srcP1, level, All );

   WALBERLA_CHECK_EQUAL( value1, value2, "values of transposes have to be equal" );

   // VTK
   VTKOutput vtkOutput( "../../output", "P1ToDGOperatorTest", storage );
   vtkOutput.add( srcDG );
   vtkOutput.add( dstP1 );
   vtkOutput.add( srcP1 );
   vtkOutput.add( dstDG );
   vtkOutput.write( level, 0 );
}

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   //checkP1ToDG1ByIntegral( 2 );
   //checkP1ToDG1ByIntegral( 3 );
   //enumerateTest();
   checkTranspose();

   return EXIT_SUCCESS;
}
