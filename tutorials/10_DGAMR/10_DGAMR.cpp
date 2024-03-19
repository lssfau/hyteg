/*
 * Copyright (c) 2022 Nils Kohl.
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

/**
 * \page 10_DGAMR Tutorial 10 - DG and AMR
 *
 * \dontinclude tutorials/10_DGAMR/10_DGAMR.cpp
 *
 * \brief In this tutorial, we solve the Poisson equation using a discontinuous Galerkin discretization and adaptive mesh
 * refinement.
 *
 * The setup is taken from <em>Elman, H. C., Silvester, D. J., & Wathen, A. J. (2014). Finite elements and fast iterative solvers:
 * with applications in incompressible fluid dynamics</em>. It is a standard test case that introduces
 * a singularity at the origin. Due to this singularity, regular global refinement does not meet the expected optimal a priori
 * error estimates. This issue is in practice usually circumvented by adaptive local refinement around the singularity.
 *
 * The analytical solution to the Poisson equation
 * \f[
    \begin{eqnarray*}
        - \Delta u = f
    \end{eqnarray*}
 * \f]
 * is given in polar coordinates by
 * \f[
    \begin{align*}
        u(r, \theta) &= r^{2/3} \sin( (2 \theta + \pi) / 3 ), \\
        f &= 0.
    \end{align*}
 * \f]
 *
 *
 * This tutorial demonstrates how adaptive mesh refinement (AMR) is applied together with the discontinuous Galerkin method.
 * For details on the discretization we refer to the literature, for example <em>Rivière, B. (2008). Discontinuous Galerkin methods
 * for solving elliptic and parabolic equations: theory and implementation</em>.
 *
 * \section T10-domainsetuprefinement Domain setup and refinement.
 *
 * We first create the PrimitiveStorage from a mesh file has has been done in previous examples.
 *
 * Mesh refinement procedures can be split into two types.
 *
 * * Conforming algorithms refine the mesh while keeping a conforming structure. That means, no hanging nodes are created. This has
 * the advantage that standard conforming finite element discretization can still be applied without any modification.
 * The disadvantage of this method is that it might produce ill-shaped elements.
 *
 * * Non-conforming algorithms may maintain the quality of the elements. However, in general those
 * methods may produce hanging nodes and the standard conforming finite element discretizations cannot be applied directly.
 * Non-conforming discretizations, on the other hand can still be applied.
 *
 * In this tutorial we cover the second case - and thus apply a discontinuous Galerkin discretization to the problem.
 * However, the first method is also applicable in HyTeG.
 *
 * The refinement is applied directly to the PrimitiveStorage and yields a tree-like structure. When a (volume-)primitive is
 * refined, child primitives are created. In 2D, faces are split into 4 face primitives, in 3D cells are split into 8 cell
 * primitives. This works as described in <em>Bey, J. (1995). Tetrahedral grid refinement</em>.
 *
 * The algorithm, however, is restricted to maintain a 2:1 balance between neighboring volume primitives. That means the
 * refinement levels of neighboring macros can only differ by 1. Before the refinement is applied, the marked primitives are
 * passed to a function that marks further primitives for refinement if necessary to keep the 2:1 balance.
 *
 * All of this happens in one method. It requires only to pass the PrimitiveID of all (local) volume primitives to be refined.
 *
 * \snippet tutorials/10_DGAMR/10_DGAMR.cpp macro refinement
 *
 * In this example we simply refine all primitives that have a vertex at the origin multiple times. The result looks like this:
 *
 * \htmlonly
    <center>
    <table>
    <tr>
    <td><img src="10_macro_0_ref.png" width="100%"/><center>initial coarse grid</center></td>
    <td><img src="10_macro_2_ref.png" width="100%"/><center>after 2 refinements</center></td>
    <td><img src="10_macro_5_ref.png" width="100%"/><center>after 5 refinements</center></td>
    </tr>
    </table>
    </center>
    \endhtmlonly
 *
 * In a second step, these coarse grids are refined uniformly.
 *
 *
 *
 * \section T10-dg Discontinuous Galerkin
 * 
 * DG functions and operators are almost used as it's done for other conforming discretizations. The main difference is that
 * the basis and degree are passed in the constructor. This also applied for the operators - the form has to be passed here, too.
 *
 * \snippet tutorials/10_DGAMR/10_DGAMR.cpp DG parameters
 *
 * \snippet tutorials/10_DGAMR/10_DGAMR.cpp DG functions
 *
 * \snippet tutorials/10_DGAMR/10_DGAMR.cpp DG operators
 *
 * The weak enforcement of Dirichlet boundary conditions is handled by a function that also takes as an argument the Form object.
 *
 * Interpolation of functions into the finite element space is currently handled by solution of the system
 * \f[
    \begin{eqnarray*}
        M u = (f, v)_\Omega
    \end{eqnarray*}
 * \f]
 * That requires the solution of a linear system involving the mass matrix \f$ M \f$:
 *
 * \snippet tutorials/10_DGAMR/10_DGAMR.cpp interpolation
 *
 * Eventually, we solve the linear system and plot the error over \f$ \frac{1}{h} \f$.
 *
 * \snippet tutorials/10_DGAMR/10_DGAMR.cpp solve
 * \snippet tutorials/10_DGAMR/10_DGAMR.cpp error
 *
 * \htmlonly
   <center>
   <table>
   <tr>
   <td><img src="10_error_conv_vs_h.png" width="100%"/><center>error vs refinement</center></td>
   <td><img src="10_error_conv_vs_dofs.png" width="100%"/><center>error vs DoFs</center></td>
   </tr>
   </table>
   </center>
   \endhtmlonly
 *
 * Standard global refinement is not sufficient to retain the optimal asymptotic accuracy. But with aggressive local refinement
 * the accuracy can be recovered.
 *
 * \htmlonly
  <center>
  <table>
  <tr>
  <td><img src="10_errors_0_2_5.png" width="100%"/><center>Plot of the error on different macro-refinement levels (0, 2, and 5). Uniform level is 3 for all three images.</center></td>
  </tr>
  </table>
  </center>
  \endhtmlonly
 *
 * \htmlonly
   <center>
   <table>
   <tr>
   <td><img src="10_sol_5_level3.png" width="100%"/><center>Plot of the computed solution and the mesh after 5 local and 3 uniform refinement steps.</center></td>
   </tr>
   </table>
   </center>
    \endhtmlonly
 */

#include <string>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

std::shared_ptr< PrimitiveStorage > buildPrimitiveStorage( uint_t refinementsAtOrigin )
{
   // Reading in a Gmsh file of an L-shaped domain.
   auto meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/lshaped.msh" );

   // Building a primitive structure from the mesh file.
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // Setting all boundary flags to 1 which defaults to Dirichlet.
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // Creating the actual, possibly distributed primitive structure.
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   /// [macro refinement]
   // Refining multiple times at the boundary.
   for ( uint_t i = 1; i <= refinementsAtOrigin; i++ )
   {
      // Here we collect all face IDs of faces that we want to refine.
      std::vector< PrimitiveID > facesToRefine;

      auto faceIDs = storage->getFaceIDs();

      // Looping over all macro-faces.
      for ( auto faceID : faceIDs )
      {
         auto face              = storage->getFace( faceID );
         auto vertexCoordinates = face->getCoordinates();

         real_t eps      = 1e-8;
         auto   atOrigin = [eps]( const Point3D& x ) { return x.norm() < eps; };

         // If any vertex of the face is near the origin, mark for refinement.
         if ( std::any_of( vertexCoordinates.begin(), vertexCoordinates.end(), atOrigin ) )
         {
            facesToRefine.push_back( faceID );
         }
      }

      // Eventually we start the refinement process. Note that possibly additional macro-faces are refined to maintain the 2:1
      // balance.
      storage->refinementHanging( facesToRefine );
   }
   /// [macro refinement]

   return storage;
}

void DGAMR( uint_t localMacroRefinements, uint_t globalMicroRefinements, std::string dbFile )
{
   using namespace dg;

   auto storage = buildPrimitiveStorage( localMacroRefinements );

   writeDomainPartitioningVTK( *storage, "./vtk", "DGAMR_Domain_MacroRef_" + std::to_string( localMacroRefinements ) );

   // r^(2/3) * sin( (2 * theta + pi) / 3 )
   auto solution = []( const Point3D& x ) {
      auto r     = x.norm();
      auto theta = std::atan2( x[1], x[0] );

      return std::pow( r, 2.0 / 3.0 ) * std::sin( ( 2 * theta + pi ) / 3 );
   };

   /// [DG parameters]
   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm = std::make_shared< DGDiffusionForm_Example >( 1, solution, solution );
   auto massForm    = std::make_shared< DGMassForm_Example >();
   /// [DG parameters]

   /// [DG functions]
   DGFunction< real_t > u( "u", storage, globalMicroRefinements, globalMicroRefinements, basis, 1 );
   DGFunction< real_t > f( "f", storage, globalMicroRefinements, globalMicroRefinements, basis, 1 );
   /// [DG functions]
   DGFunction< real_t > sol( "sol", storage, globalMicroRefinements, globalMicroRefinements, basis, 1 );
   DGFunction< real_t > tmp( "tmp", storage, globalMicroRefinements, globalMicroRefinements, basis, 1 );
   DGFunction< real_t > Merr( "Merr", storage, globalMicroRefinements, globalMicroRefinements, basis, 1 );
   DGFunction< real_t > err( "err", storage, globalMicroRefinements, globalMicroRefinements, basis, 1 );

   DGFunction< idx_t > numerator( "numerator", storage, globalMicroRefinements, globalMicroRefinements, basis, 1 );
   numerator.enumerate( globalMicroRefinements );

   /// [DG operators]
   DGOperator A( storage, globalMicroRefinements, globalMicroRefinements, laplaceForm );
   DGOperator M( storage, globalMicroRefinements, globalMicroRefinements, massForm );
   /// [DG operators]

   // Assemble RHS. RHS == 0.
   f.applyDirichletBoundaryConditions( laplaceForm, globalMicroRefinements );

   /// [interpolation]
   tmp.evaluateLinearFunctional( solution, globalMicroRefinements );

   PETScCGSolver< DGOperator > solverM( storage, globalMicroRefinements, numerator );

   solverM.solve( M, sol, tmp, globalMicroRefinements );
   /// [interpolation]

   // Solve system.
   /// [solve]
   PETScCGSolver< DGOperator > solverA( storage, globalMicroRefinements, numerator, 1e-12, 1e-12, 10000 );
   solverA.solve( A, u, f, globalMicroRefinements );
   /// [solve]

   /// [error]
   err.assign( { 1.0, -1.0 }, { u, sol }, globalMicroRefinements );

   M.apply( err, Merr, globalMicroRefinements, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, globalMicroRefinements, Inner ) );
   /// [error]

   auto numDofs = u.getNumberOfGlobalDoFs( globalMicroRefinements );

   WALBERLA_LOG_INFO_ON_ROOT( "Number of DoFs: " << numDofs );
   WALBERLA_LOG_INFO_ON_ROOT( "Error (macro-ref: " << localMacroRefinements << " | level: " << globalMicroRefinements
                                                   << "): " << discrL2 );

   VTKOutput vtk( "./vtk/", "DGAMR_Functions_MacroRef_" + std::to_string( localMacroRefinements ), storage );
   vtk.add( u );
   vtk.add( sol );
   vtk.add( err );
   vtk.write( globalMicroRefinements );

   FixedSizeSQLDB db( dbFile );

   db.setVariableEntry( "local_refinements", localMacroRefinements );
   db.setVariableEntry( "global_refinements", globalMicroRefinements );
   db.setVariableEntry( "num_dofs", numDofs );
   db.setVariableEntry( "error_l2", discrL2 );

   db.writeRowOnRoot();
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   std::string dbFile = "dgamr.db";

   DGAMR( 0, 1, dbFile );
   DGAMR( 0, 2, dbFile );
   DGAMR( 0, 3, dbFile );
   DGAMR( 0, 4, dbFile );
   DGAMR( 0, 5, dbFile );
   DGAMR( 0, 6, dbFile );
   DGAMR( 0, 7, dbFile );

   DGAMR( 2, 1, dbFile );
   DGAMR( 2, 2, dbFile );
   DGAMR( 2, 3, dbFile );
   DGAMR( 2, 4, dbFile );
   DGAMR( 2, 5, dbFile );
   DGAMR( 2, 6, dbFile );

   DGAMR( 5, 1, dbFile );
   DGAMR( 5, 2, dbFile );
   DGAMR( 5, 3, dbFile );
   DGAMR( 5, 4, dbFile );
   DGAMR( 5, 5, dbFile );
   DGAMR( 5, 6, dbFile );
}
