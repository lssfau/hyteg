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

/**
 * \page 09_BoundaryConditions Tutorial describing use of multiple boundary conditions.
 *
 * \dontinclude tutorials/09_BoundaryConditions/09_BoundaryConditions.cpp
 *
 * \brief In this tutorial, we are going to demonstrate how one can set different types of boundary conditions
 * on different parts of the boundary and/or for different components of the problem.
 *
 * \section Task
 * Assume that we want to solve a convection problem for a strongly viscous fluid. We assume that convection
 * is driven by thermally induced density differences and that a Boussinesq approximation is employed. This
 * will leads us, in a simple case, to a PDE system similar to the following one:
 * \f[
   \begin{eqnarray*}
   - \Delta \underline{u} + \nabla p &=& \underline{F} \\[2ex]
   \nabla\cdot \underline{u} &=& 0 \\[1.5ex]
   \frac{\partial T}{\partial t} + \underline{u}\cdot\nabla T - \Delta \underline T &=& \underline{G}
   \end{eqnarray*}
 * \f]
 * Here \f$\underline{u}(\underline{x},t),p(\underline{x},t),T(\underline{x},t)\f$ are our unknown functions
 * representing velocity, pressure and temperature. We pose the problem on a rectangular domain
 * \f$\Omega = (0,1.5)\times(0,1)\f$ and need to add boundary conditions for \f$\underline{u}\f$ and \f$T\f$
 * at the boundaries of our rectangle.
 * 
 * For demonstration purposes we use the following boundary conditions:
 * <table border="0" width="90%"><tr>
 * <td><img src="09_BoundaryConditions_Velocity.png" width="95%" /></td>
 * <td><img src="09_BoundaryConditions_Temperature.png" width="95%" /></td>
 * </tr></table>
 * Expressed mathematically we require for the velocity
 * \f[
   \begin{array}{lll}
   \mbox{bottom edge:} &
   \underline{u}(x,y,t) = (0,0)^\top & \forall (x,0), x\in(0,1.5) \\[2ex]
   \mbox{top edge:} &
   \underline{u}(x,y,t) = (2,0)^\top & \forall (x,1), x\in(0,1.5) \\[2ex]
   \mbox{left and right edge:} & \displaystyle
   \underline{u}\cdot\hat{n} = 0 \mbox{ and }
   \frac{\partial u_y}{\partial \hat{n}} (x,y,t) = 0 &
   \forall (x,y), x\in\{0,1.5\}, y\in(0,1)
   \end{array}
 * \f]
 * where \f$\hat{n}\f$ represents an outward normal, and similarly for the temperature
 * \f[
   \begin{array}{lll}
   \mbox{bottom edge:} &
   T(x,y,t) = 1 & \forall (x,0), x\in(0,1.5) \\[2ex]
   \mbox{top edge:} &
   T(x,y,t) = 0 & \forall (x,1), x\in(0,1.5) \\[2ex]
   \mbox{left and right edge:} &
   \displaystyle\frac{\partial T}{\partial \hat{n}} (x,y,t) = 0 &
   \forall (x,y), x\in\{0,1.5\}, y\in(0,1)
   \end{array}
 * \f]
 * which both must hold for all times \f$t\in(t_0,t_1)\f$.
 *
 * \section Implementation
 * 
 * We will now take a look at how we can set these boundary conditions in HyTeG.
 *
 * \subsection Step1 Step #1: Mark parts of the Boundary
 * 
 * The first step we need to perform is mark the four different parts of the boundary of our problem domain,
 * i.e. the four edges of our rectangle, by setting corresponding **MeshBoundaryFlags**. These are simple
 * integer values that will allows us to differentiate between the parts when setting the boundary conditions.
 * We are going to use the following values:
 *
 * <img src="09_BoundaryConditions_Flags.png" width="45%" /></td>
 * 
 * The flags will be set on the \link hyteg::SetupPrimitiveStorage `SetupPrimitiveStorage`\endlink object we
 * generate from our mesh using one of the different `setMeshBoundaryFlags` methods. In our case
 * \link hyteg::SetupPrimitiveStorage::setMeshBoundaryFlagsByCentroidLocation() `setMeshBoundaryFlagsByCentroidLocation()`\endlink
 * is the appropriate one.
 * We need to supply to it a callback that returns true or false, depending on whether the centroid of a macro primitive
 * belongs to the boundary part we want to flag, or not.
 * Note that we might also use 
 * \link hyteg::SetupPrimitiveStorage::setMeshBoundaryFlagsByVertexLocation() `setMeshBoundaryFlagsByVertexLocation()`\endlink.
 * However, this method requires a little more consideration w.r.t. the callbacks as the flags for higher-dimensional
 * primitives are derived from those of their associated vertices.
 * 
 * \snippet tutorials/09_BoundaryConditions/09_BoundaryConditions.cpp Flag_Boundaries_Parts
 *
 * \subsection Step2 Step #2: Create BoundaryCondition objects
 * 
 * The next step is to create one object of type \link hyteg::BoundaryCondition `BoundaryCondition` \endlink
 * for each unknown function we have in our
 * problem. In our case we will generate two, one for the velocity and one for the temperature. Note that
 * we do not explicitely create one for pressure, a detail to which we will come back below.
 *
 * Once we created an object, we can make it aware of what type of boundary condition (Dirichlet, Neumann, ...)
 * is required for this type of function on which part of the boundary. Here the flags (integer values) we set
 * in step 1 come back into play.
 * The methods we use here return a \link hyteg::BoundaryUID `BoundaryUID` \endlink. One
 * way to set the actual boundary values uses this information, see Step #4. If we do not need the UIDs we
 * can just ignore the return value. In our example we do this for the Neumann and free-slip parts.
 * 
 * Note that the symbolic names we supply as arguments are stored internally in the object and must be unique,
 * but currently find no further usage.
 * 
 * \snippet tutorials/09_BoundaryConditions/09_BoundaryConditions.cpp BC_Objects
 * 
 * \subsection Step3 Step #3: Create Functions
 * 
 * Now we can instantiate function objects for our unknown functions. We are going to use a P<sub>2</sub> space for temperature
 * and velocity and a P<sub>1</sub> space for pressure. Note that in HyTeG each function knowns about its boundary condition.
 * 
 * We can provide this information to its constructor by passing a corresponding `BoundaryCondition`
 * object. This is what we are going to do in the case of `temperature` and `velocity` below. Alternatively we can use the
 * setter method `setBoundaryCondition()` on the function object once it was created.
 *
 * \note If we do not pass a `BoundaryCondition` object to the constructor, a default one will be used. The constructor uses
 * the method \link hyteg::BoundaryCondition::create0123BC() `create0123BC()`\endlink for this. The latter will set the
 * following types
 * <center>
 * <table>
 * <tr><td align="center">mesh boundary flag</td><td align="center">BC</td></tr>
 * <tr><td align="right">    1</td><td align="left">Dirichlet</td></tr>
 * <tr><td align="right">    2</td><td align="left">Neumann</td></tr>
 * <tr><td align="right">    3</td><td align="left">Free-slip</td></tr>
 * <tr><td align="right">other</td><td align="left">Inner</td></tr>
 * </table>
 * </center>
 * This is the reason, why in previous tutorials, where only one type
 * of boundary condition was present, mesh boundary flags were set using 
 * \link hyteg::SetupPrimitiveStorage::setMeshBoundaryFlagsOnBoundary() `setMeshBoundaryFlagsOnBoundary()`\endlink and no
 * `BoundaryCondition` object was explicitely created.
 * As an example consider the following call from tutorial #8.
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp create_storage neumann
 * It marks all parts of the boundary by assigning their vertices a flag
 * value of 2, and all other vertices by flagging them with 0. In combination with the object returned by `create0123BC()`
 * this corresponds to Neumann boundary conditions everywhere on the boundary.
 *
 * For `pressure` the problem specifies no boundary conditions at all. Thus, we set all flags to be 'Inner' by constructing
 * a `BoundaryCondition` object using  \link hyteg::BoundaryCondition::createAllInnerBC() `createAllInnerBC()`\endlink.
 *
 * Normally you would not generate the velocity and pressure as two different functions, but instead use a
 * composite function, i.e. in our case a `P2P1TaylorHoodFunction`. The latter knows that its second scalar
 * function has no boundary conditions. Thus, the constructor will use the passed `BoundaryCondition` object
 * only for the velocity component. For the pressure component an object will be generated automatically via
 * \link hyteg::BoundaryCondition::createAllInnerBC() `createAllInnerBC()`\endlink.
 * 
 * \snippet tutorials/09_BoundaryConditions/09_BoundaryConditions.cpp Function_Creation
 * 
 * \subsection Step4 Step #4: Setting Boundary and Initial Conditions
 *
 * Setting values on our function is done by using the `%interpolate()` method. There are two possible ways
 * to specify which degrees of freedom we want to change:
 *
 * 1. Each degree of freedom has a type attached. Thus, we can specify a flag of type `hyteg::DoFType` to
 *    select certain degrees of freedom. The possible values are `All`, `Boundary`, `Inner`, `DirichletBoundary`,
 *    `NeumannBoundary`, `FreeslipBoundary`. These can be combine by `|`.  
 *
 *    Note that by calling `bcTemperature.createNeumannBC( "neumannWalls", { 3, 4 } )` in step #2 and passing
 *    `bcTemperature` to the constructor of our temperature function in step #3, we have set the DoFType of
 *    the temperature degrees of freedom on the left and right boundary to be `NeumannBoundary`.
 * 2. The other way to select degrees of freedom is to use a corresponding `BoundaryUID`. By this we can
 *    directly address degrees of freedom of a function belonging to a specific boundary condition.
 *
 * Let us set the Dirichlet values for temperature now. Normally, if the value we want to set is a constant,
 * we can directly use it without the need for an `std::function`. However, our problem features two different
 * Dirichlet values on top and bottom. As there is no way to express this via the DoFType flag, we have to write
 * a short lambda expression:
 *
 * \snippet tutorials/09_BoundaryConditions/09_BoundaryConditions.cpp Setting_BC_Temp
 *
 * When we set boundary values for velocity we must to take into account that this is a vector-valued function.
 * Thus, we need to either specify one constant or one std::function per component.
 *
 * \snippet tutorials/09_BoundaryConditions/09_BoundaryConditions.cpp Setting_BC_Vel
 *
 * As a final aspect let us demonstrate how to combine different DoFType values to set the initial values for
 * the temperature. Note that the initiation function
 *
 * \f[
 (x,y) \mapsto (1-y)+\frac{1}{100}\cos\Big(\frac{\pi x}{1.5}\Big)\sin\Big(\pi y\Big) 
 * \f]
 * is consistent with the Dirichlet boundary conditions so we could also use `All` to combine setting boundary
 * and initial conditions for time \f$t=t_0\f$.
 * \snippet tutorials/09_BoundaryConditions/09_BoundaryConditions.cpp Setting_ICs
 *
 */

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::Environment env( argc,
                              argv ); // TODO: this needs to be added for MPI builds I think (did not build on my machine ;))
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t nx       = 3;
   const uint_t ny       = 2;
   const uint_t minLevel = 2;
   const uint_t maxLevel = 3;

   /// [Flag_Boundaries_Parts]
   // generate mesh and setup storage
   auto meshInfo     = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1.5, 1} ), MeshInfo::CROSS, nx, ny );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // The following statement is not strictly necessary. It sets the mesh boundary flags to a well-defined state,
   // overwriting whatever came in via the MeshInfo object.
   //
   // In our case the MeshInfo::meshRectangle() function creates a MeshInfo object where vertices and edges on
   // the boundary are marked as 1 and vertices, edges and faces in the interior as 0, anyway.
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // callbacks for marking parts of the boundary
   real_t eps         = 1e-3;
   auto   topBoundary = [eps]( const Point3D& x ) { return std::abs( x[1] - real_c( 1 ) ) < eps; };

   auto bottomBoundary = [eps]( const Point3D& x ) { return x[1] < eps; };

   auto leftBoundary = [eps]( const Point3D& x ) {
      return x[0] < eps && !( std::abs( x[1] - real_c( 1 ) ) < eps ) && !( x[1] < eps );
   };

   auto rightBoundary = [eps]( const Point3D& x ) {
      return std::abs( x[0] - real_c( 1.5 ) ) < eps && !( std::abs( x[1] - real_c( 1 ) ) < eps ) && !( x[1] < eps );
   };

   // assigning mesh boundary flags to different parts of the boundary
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 1, topBoundary );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 2, bottomBoundary );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 3, leftBoundary );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 4, rightBoundary );

   // now create the actual storage
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );
   /// [Flag_Boundaries_Parts]

   /// [BC_Objects]
   BoundaryCondition bcVelocity;
   BoundaryCondition bcTemperature;

   BoundaryUID idTopWallVel = bcVelocity.createDirichletBC( "topWall", {1} );
   BoundaryUID idBotWallVel = bcVelocity.createDirichletBC( "botWall", {2} );
   bcVelocity.createFreeslipBC( "sideWalls", {3, 4} );

   BoundaryUID idTopWallTemp = bcTemperature.createDirichletBC( "dirichletWallTop", {1} );
   BoundaryUID idBotWallTemp = bcTemperature.createDirichletBC( "dirichletWallBot", {2} );
   bcTemperature.createNeumannBC( "neumannWalls", {3, 4} );
   /// [BC_Objects]

   /// [Function_Creation]
   P2Function< real_t > temperature( "temperature", storage, minLevel, maxLevel, bcTemperature );

   // for demonstration purposes create u and p separately ...
   P2VectorFunction< real_t > velocity( "velocity", storage, minLevel, maxLevel, bcVelocity );
   P1Function< real_t >       pressure( "pressure", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() );

   // ... normally you want to use a combined function
   P2P1TaylorHoodFunction< real_t > uAndp( "stokes function", storage, minLevel, maxLevel, bcVelocity );
   /// [Function_Creation]

   /// [Setting_BC_Temp]
   // fixed temperatures on top and bottom edge
   real_t topT = real_c( 0.0 );
   real_t botT = real_c( 1.0 );

   // as the values are a constant we do not need to provide a callback function to the interpolation method
   temperature.interpolate( botT, maxLevel, idBotWallTemp );
   temperature.interpolate( topT, maxLevel, idTopWallTemp );
   /// [Setting_BC_Temp]

   /// [Setting_BC_Vel]

   // for a vector-valued function we need one std::function or constant
   // per component (note that we cannot mix the types, though)

   // set velocity field on bottom wall to no-slip
   velocity.interpolate( {real_c( 0 ), real_c( 0 )}, maxLevel, idBotWallVel );

   // set velocity field on top wall to no-outflow + non-zero tangential component
   velocity.interpolate( {real_c( 2 ), real_c( 0 )}, maxLevel, idTopWallVel );

   // to set the boundary values on our TaylorHood funcion, we extract its velocity component
   uAndp.uvw().interpolate( {real_c( 0 ), real_c( 0 )}, maxLevel, idBotWallVel );
   uAndp.uvw().interpolate( {real_c( 2 ), real_c( 0 )}, maxLevel, idTopWallVel );

   /// [Setting_BC_Vel]

   /// [Setting_ICs]
   auto initialTemperature = []( const Point3D& x ) {
      return ( 1.0 - x[1] ) + 0.01 * std::cos( pi * x[0] / real_c( 1.5 ) ) * std::sin( pi * x[1] );
   };
   temperature.interpolate( initialTemperature, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );
   /// [Setting_ICs]

   // export functions to visually check results
   bool beVerbose = true;
   if ( beVerbose )
   {
      std::string fPath = ".";
      std::string fName = "boundaryValues";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( velocity );
      vtkOutput.add( temperature );
      vtkOutput.add( pressure );
      vtkOutput.add( uAndp );
      vtkOutput.write( maxLevel );
   }
}
