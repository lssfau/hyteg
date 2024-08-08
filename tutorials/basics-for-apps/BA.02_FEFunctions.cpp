/*
 * Copyright (c) 2024 Marcus Mohr.
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
#include "core/mpi/Environment.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using namespace hyteg;

/**
 * \page BA.02_FEFunctions Tutorial BA.02 - Setting up FE Functions
 *
 * \dontinclude tutorials/basics-for-apps/BA.02_FEFunctions.cpp
 *
 * \brief In this tutorial we will demonstrate how to set up a C++ object
 * that represents a Finite Element function.
 *
 * \section BA02-intro Introduction
 *
 * A large variety of function spaces exists for performing Finite Element
 * simulations. HyTeG supports several of these. In this tutorial we are
 * going to learn how to generate a Finite Element function and how to
 * set its values.
 * 
 * \section BA02-mesh Preparing the PrimitiveStorage
 * 
 * In order to generate a function, we must first prepare a domain
 * representation in our application. As explained in \ref BA.01_PrimitiveStorage
 * this finally leads us to a PrimitiveStorage. The basic workflow is as before.
 * However, this time an additional step is involved:
 *
 * \snippet{trimleft} this FlaggingDoFs
 * 
 * In an actual simulation one often wants to initialise the degrees of freedom (DoFs)
 * of a function differently depending on their type. A typical example would be
 * a PDE problem with Dirichlet boundary conditions for which we create an FE function
 * object that will later represent the discrete weak solution of the problem. In
 * this case we might want to intialise all DoFs on the Dirichlet boundary to the
 * values given by the boundary condition function and those in the interior of the domain to zero.
 *
 * As a first step to be able to do this we call SetupPrimitiveStorage::setMeshBoundaryFlagsOnBoundary.
 * This member function will mark all primitives belonging to the domain boundary with the value of the
 * first argument (1 in our example) and those in the interior with the second (0 in our example).
 * The final third argument indicates that the primitives with highest geometric dimension, e.g.
 * faces in a 2D mesh, can automatically be considered to belong to the interior of the domain.
 * 
 * A second small difference is, that this time we do not use a PrimitiveStorage object directly,
 * but instead a shared pointer to it.
 *
 * \section BA02-P1Function Using a P1 Function to represent a scalar field
 * 
 * Now let us assume that we need to represent on our domain a scalar field, say pressure as a
 * physical quantity. One of the simplest kinds of Finite Element functions is the class of
 * continuous Lagrange elements of order 1, so continuous piecewise linear polynomials.
 * In HyTeG this is represented by the hyteg::P1Function class. (Note that the latter is
 * actually only an alias for the hyteg::VertexDoFFunction class.)
 * 
 * Setting up a corresponding object is a one-liner:
 * 
 * \snippet{trimleft} this SimpleP1Function
 *
 * The arguments we pass have the following meaning:
 * - The first is the name of the function. This will e.g. be used when we export data for
 *   visualisation. It is also used by HyTeG internally at some places to distinguish or
 *   identify functions. Thus, the names you choose should be globally unique within your
 *   application.
 * - A Finite Element function is inevitably tied to the underlying spatial discretisation.
 *   Thus, the second argument is our smart pointer to the PrimitiveStorage object.
 * - The Hybrid Hierachical Grids paradigm employed by HyTeG constructs a sequence of meshes
 *   starting from a coarse base mesh. This provides a nested hierarchy of meshes, which can
 *   e.g. be used by geometric multigrid methods. Hence, the third and fourth argument we
 *   provide specify the range of levels from coarsest to finest level on which the function "lives".
 *   It is completely fine for these to be the same, like e.g. the finest mesh level on which
 *   you intend to solve the PDE problem.
 *
 * Also noteworthy is that all Finite Element functions in HyTeG are templated with the datatype
 * used for the degrees of freedom. In our case we use walberla::real_t, which is HyTeG's internal
 * floating point data type (inherited from the underlying waLBerla library). It defaults to double.
 * But this can be changed to float when running CMake by setting WALBERLA_DOUBLE_ACCURACY=no.
 *
 * After creation all degrees of freedom of our P1Function are zero on all levels. Let us say, that
 * we want to change this and make the P1Function interpolate the mathematical function given as
 *
\f[ (x,y) \mapsto 2x^2 + \frac{1}{2} y^3 \f]
 *
 * In order to achieve this we need to prepare a <a href="https://en.cppreference.com/w/cpp/utility/functional/function">std::function</a> object.
 * The latter is a C++ class that can represent various different callables. For our purposes using a lambda function is
 * the most convenient:
 *
 * \snippet{trimleft} this Interpolation1
 *
 * Note that our lambda function needs to accept coordinates for evaluation in the form of hyteg::Point3D,
 * although we are in 2D. walberla::real_c takes care of converting the literals to the real_t datatype.
 * This is not necessary, if you directly use literals of the corresponding type or can live with implicit
 * type conversion.
 *
 * Now that we have created an "expression" that represents our mathematical function, we can use it to
 * set the DoF values of our Finite Element function. For this we call the interpolate member function
 * and pass the expression. The next argument gives the level on which we want to set the DoF values.
 * The final one tells the function to set all values independent of their hyteg::DoFType.
 * 
 * \snippet{trimleft} this Interpolation2
 * 
 * If we want to set the DoFs to a constant value, we can do so directly, without the need for a std::function.
 * Let us set all DoFs on the coarsest level on which the pressure function was defined and which belong to
 * the domain interior to a value of 2.35. We can achieve this via:
 * 
 * \snippet{trimleft} this Interpolation3
 *
 * \section BA02-P2VectorFunction Representing a vector field
 * 
 * In Finite Element simulations one often needs to represent not only scalar, but also vector fields.
 * Think for example on the displacement field in an elasticity problem, the velocity field in
 * computational fluid dynamics, or the magnetic field in an electrodynamics problem.
 * 
 * Let us take a look here at the straightforward approach of representing a vector field by simply using
 * scalar Finite Element functions for its components. In HyTeG the corresponding function classes are
 * children of hyteg::CSFVectorFunction.
 *
 * Let us extend our example by assuming that we also have a velocity field in our problem and want to
 * represent this with quadratic Lagrange elements. For this we can use the hyteg::P2VectorFunction class.
 * Setting up a corresponding object works as in the scalar case:
 *
 * \snippet{trimleft} this P2VectorFunction
 *
 * Note that the constructor will automatically determine the number of components from the PrimitiveStorage
 * we pass. If this contains macro-cells the function will have three components, if the highest dimensional
 * primitive is a macro-face, it will have two.
 *
 * As we are basically using a continuous Lagrange element, we can set the DoFs again by interpolation. So,
 * let us define two expressions, one for the x-component and one for the y-component of velocity. Then we
 * use these to set the DoFs which belong to the Dirichlet boundary of the domain:
 *
 * \snippet{trimleft} this Interpolation4
 *
 * We can use constant values directly. However, currently it is not possible to mix std::functions and
 * constant values in one call to interpolate. As an alternative we can set the values of the different
 * component functions separately by requesting a reference from the P2VelocityFunction object. As an
 * example let us change the DoF values of the y-component (index 1) to be 3.0 on the boundary:
 *
 * \snippet{trimleft} this Interpolation5
 *
 * This concludes our tutorial on setting up Finite Element functions. Take a look at \ref BA.03_ExportForVisualisation
 * to learn how you can export them for visualisation.
 *
 * \section BA02-Code Complete Program
 * \include tutorials/basics-for-apps/BA.02_FEFunctions.cpp
 * 
 */

int main( int argc, char** argv )
{
   walberla::mpi::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   /// [FlaggingDoFs]
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "2D/bfs_12el.msh" ) );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // new step
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   /// [FlaggingDoFs]

   const uint_t minLevel = 2;
   const uint_t maxLevel = 4;

   /// [SimpleP1Function]
   P1Function< real_t > pressure( "Pressure", storage, minLevel, maxLevel );
   /// [SimpleP1Function]

   /// [Interpolation1]
   std::function< real_t( const hyteg::Point3D& ) > pressureExpression = []( const hyteg::Point3D& x ) {
      return real_c( 2 ) * x[0] * x[0] + real_c( 0.5 ) * x[1] * x[1] * x[1];
   };
   /// [Interpolation1]

   /// [Interpolation2]
   uint_t level = maxLevel;
   pressure.interpolate( pressureExpression, level, All );
   /// [Interpolation2]

   /// [Interpolation3]
   pressure.interpolate( real_c( 2.35 ), minLevel, Inner );
   /// [Interpolation3]

   /// [P2VectorFunction]
   P2VectorFunction< real_t > velocity( "Velocity Field", storage, minLevel, maxLevel );
   /// [P2VectorFunction]

   /// [Interpolation4]
   std::function< real_t( const hyteg::Point3D& ) > xExpr = []( const hyteg::Point3D& x ) {
      return real_c( 2 ) * x[0] * x[0] - x[1];
   };

   std::function< real_t( const hyteg::Point3D& ) > yExpr = []( const hyteg::Point3D& x ) {
      return real_c( -1.5 ) * x[0] + real_c( 2.0 ) * std::sin( x[1] );
   };

   velocity.interpolate( { xExpr, yExpr }, level, DirichletBoundary );
   /// [Interpolation4]

   /// [Interpolation5]
   velocity[1].interpolate( real_c( 3 ), level, DirichletBoundary );
   /// [Interpolation5]
}
