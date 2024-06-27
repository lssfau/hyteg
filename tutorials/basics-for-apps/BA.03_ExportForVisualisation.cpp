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

#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using namespace hyteg;

/**
 * \page BA.03_ExportForVisualisation Tutorial BA.03 - Exporting FE Functions for Visualisation
 *
 * \dontinclude tutorials/basics-for-apps/BA.03_ExportForVisualisation.cpp
 *
 * \brief In this tutorial we will demonstrate how one can export FE functions from a simulation
 * to visualise them.
 *
 * \section BA03-preps Preparations
 *
 * In order to demonstrate how to export data for later visualisation, we need to setup
 * some finite element functions. In order to not use completely virtual data we mimick
 * solving a simple Poiseuille flow problem, i.e. laminar flow through a channel. As
 * domain we choose a rectangle \f$ \Omega = (0,4) \times (-1,1)\f$ and use one of HyTeG's
 * internal mesh generators to discretise the domain with a mesh 
 *
 * \snippet{trimleft} this PrepareDomain
 *
 * Our parameter choices for MeshInfo::meshRectangle result in the following macro mesh
 *
 * <img src="BA.03_Channel_Coarse_Mesh.png" width="65%" />
 *
 * We do not truly solve the problem but instead describe an linear pressure drop function
 * given by
 *
\f[ p : (x,y) \mapsto 1 - \frac{x}{4} \f]
 *
 * and a parabolic velocity profile
 *
\f[ u : (x,y) \mapsto \left( 1 - y^2, 0 \right)^\top \f]
 *
 * Below you see the code for generating corresponding std::function objects
 *
 * \snippet{trimleft} this MimickSolving
 *
 * Now, in the way we have seen in the previous tutorials, we setup two finite element functions,
 * a \f$P_1\f$ function to represent pressure and a \f$P_2\f$ vector function to represent the
 * velocity and interpolate our expressions.
 *
 * \snippet{trimleft} this FunctionSetup
 *
 * Note: When actually simulating the problem we would better use a hyteg::P2P1TaylorHoodFunction
 *       instead.
 * 
 * Now we are prepared to export our functions.
 *
 * \section BA03-vtkExport Exporting Data with VTKOutput
 *
 * HyTeG natively, i.e. without linking against external libraries, supports exporting finite
 * element functions and the associated mesh information in the form of <a href="https://vtk.org/">vtu-files</a>.
 * These can then be used for visualisation e.g. with <a href="https://www.paraview.org/">ParaView</a>.
 *
 * In order to do so, we first need to generate an object of type hyteg::VTKOutput
 *
 * \snippet{trimleft} this VTKOutputObject
 *
 * The arguments to the constructor call are
 *
 * - a path to the output directory where files should be placed
 * - a basename for the output files
 * - a shared pointer to a PrimitiveStorage object
 *
 * Once we have generated the object, we can register those functions with it, which we want to
 * include into our data export
 *
 *  \snippet{trimleft} this VTKOutputAdd
 *
 * Then, when we want to perform the data export, we simply call write():
 *
 *  \snippet{trimleft} this VTKOutputWrite
 *
 * This will export data for all functions that are currently registered with the object.
 * As our functions can live on multiple refinement levels we need to pass the information for
 * which level we want to export. For this to work, all functions that we added need to be defined
 * on the requested level.
 *
 * Now let us compile this tutorial and run it with two MPI processes
 *
 * \verbatim
mpirun -n 2 ./BA.03_ExportForVisualisation
\endverbatim
 *
 * As a result we obtain in the current working directory (since we specified ".") the following
 * two files:
 *
 * - channelFlow_VertexDoF_level3_ts0.vtu
 * - channelFlow_P2_level3_ts0.vtu
 *
 * The first one contains the data of our \f$P_1\f$ pressure function and the second that of our
 * \f$P_2\f$ velocity function.
 *
 * HyTeG will generate one file for each family of functions that can be represented, from the
 * perspective of VTK, on the same mesh. Thus, if we were to export two \f$P_1\f$ functions, they
 * would be written into the same file. However, as the \f$P_2\f$ function has more degrees of freedom,
 * whose (virtual) coordinates we need to export for visualisation, we need to put it into another
 * file. However, there will always be only one file per family of functions, independent of the
 * number of MPI processes.
 *
 * As one can see, to the basename we provided, which was **channelFlow**, a tag for the family of
 * functions as well as the refinement level and a marker for the time-step, a.k.a. data export
 * episode was appended.
 *
 * Since we exported data for `maxLevel = 3` these are for the thrice refined macro mesh:
 *
 * <img src="BA.03_Channel_Coarse+Fine_Mesh.png" width="65%" />
 *
 * If we want to export data for multiple time-steps, we pass a corresponding integer value
 * to the write member function:
 *
 * \code{.cpp}
vtkWriter::write( maxLevel, stepCounter );
\endcode
 *
 * For further details, see the documentation of hyteg::VTKOutput.
 *
 * Now that we have exported the data, we can visualise them e.g. with ParaView:
 * <table border="0">
 * <tr><td>
 * <img src="BA.03_Velocity.png" width="65%" />
 * </td></tr>
 * <tr><td>
 * Velocity Field
 * </td></tr>
 * <tr><td>
 * <img src="BA.03_Pressure.png" width="65%" />
 * </td></tr>
 * <tr><td>
 * Pressure
 * </td></tr>
 * </table>
 *
 * Note: It is possible to deregister functions by using hyteg::VTKOutput::remove().
 * They will then not be included into future data export episodes.
 *
 * \section BA03-adiosExport Exporting Data with ADIOS2
 *
 * As an alternative way to output data HyTeG sports an interface to
 * <a href="https://adios2.readthedocs.io/en/latest/">ADIOS 2: The Adaptable Input/Output System version 2</a>.
 * This library whose development was funded by the Exascale Computing Project (ECP) of the U.S. Department of
 * Energy constitutes a unified high-performance framework for extreme scale I/O.
 *
 * In order to make use of ADIOS2 it, naturally, must be installed on your target system. Additionally while
 * running CMake to configure HyTeG you need to specify `-DHYTEG_BUILD_WITH_ADIOS2=yes`.
 *
 * The class that your app will interface with is hyteg::AdiosWriter. Conceptually this works very similar to
 * `VTKOutput` described above. Thus, we have to
 * - create an AdiosWriter object
 * - register functions with it
 * - call `write()` to export these for a given refinement level
 *
 * \snippet{trimleft} this AdiosWriter
 *
 * Running our tutorial with Adios support enabled we get two output "files"
 *
 * - channelFlow-P1_level3.bp
 * - channelFlow-P2_level3.bp
 *
 * These are actually directories containing multiple files in which our data and some meta-information is
 * stored. Adios provides a tool `bpls` to inspect the contents of such a "BP file". It is quite powerful
 * offering many different option. A minimal query with `bpls -a channelFlow-P2_level3.bp` gives the following
 * output
 * \verbatim
  int64_t   NumberOfElements  {2}
  int64_t   NumberOfVertices  {2}
  string    Software          attr
  double    TIME              scalar
  double    Velocity Field    [2]*{918, 2}
  int64_t   connectivity      [2]*{384, 7}
  uint32_t  types             scalar
  double    vertices          [2]*{918, 3}
  string    vtk.xml           attr
\endverbatim
 *
 * Note that, as with the VTKOutput, we get two different files for the \f$P_1\f$ and  \f$P_2\f$ family
 * of functions. This is again for the sake of VTK, as they have different numbers and locations of DoFs.
 * However, using the AdiosWriter we only get one channelFlow-P1_level3.bp file independent of how many
 * export episodes, i.e. time-steps or iterations, we perform. Hence, the mesh information is stored only
 * once.
 *
 * A further advantage is that Adios can be configured, e.g. via passing a parameter file, to tune various
 * aspects of the IO. One can, e.g. configure it such that all MPI processes on a node of a distributed
 * system agglomerate their data onto a single process on that node, which then performs the I/O. For full
 * details please take a look at the <a href="https://adios2.readthedocs.io/en/latest/">documentation</a>.
 *
 * Note: In order to import the "BP files" into ParaView, select the ADIOS2VTXReader.
 *
 * \section BA03-Code Complete Program
 * \include tutorials/basics-for-apps/BA.03_ExportForVisualisation.cpp
 *
 */

int main( int argc, char** argv )
{
   walberla::mpi::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   /// [PrepareDomain]
   MeshInfo              meshInfo = MeshInfo::meshRectangle( Point2D( 0.0, -1.0 ), Point2D( 4.0, 1.0 ), MeshInfo::DIAMOND, 3, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage );
   /// [PrepareDomain]

   // specialised function to export data of the macro mesh
   hyteg::writeDomainPartitioningVTK( storage, ".", "channel_macro_mesh" );

   const uint_t minLevel = 3;
   const uint_t maxLevel = 3;

   /// [MimickSolving]
   std::function< real_t( const hyteg::Point3D& ) > pressureExpr = []( const hyteg::Point3D& x ) {
      return real_c( 1.0 ) - real_c( 0.25 ) * x[0];
   };

   std::function< real_t( const hyteg::Point3D& ) > xVelocityExpr = []( const hyteg::Point3D& x ) {
      return real_c( 1.0 ) - x[1] * x[1];
   };
   /// [MimickSolving]

   /// [FunctionSetup]
   P1Function< real_t >       pressure( "Pressure", storage, minLevel, maxLevel );
   P2VectorFunction< real_t > velocity( "Velocity Field", storage, minLevel, maxLevel );

   pressure.interpolate( pressureExpr, maxLevel );
   velocity[0].interpolate( xVelocityExpr, maxLevel );
   velocity[1].interpolate( real_c( 0 ), maxLevel );
   /// [FunctionSetup]

   /// [VTKOutputObject]
   VTKOutput vtkWriter( ".", "channelFlow", storage );
   /// [VTKOutputObject]

   /// [VTKOutputAdd]
   vtkWriter.add( pressure );
   vtkWriter.add( velocity );
   /// [VTKOutputAdd]

   /// [VTKOutputWrite]
   vtkWriter.write( maxLevel );
   /// [VTKOutputWrite]

#ifdef HYTEG_BUILD_WITH_ADIOS2
   /// [AdiosWriter]
   AdiosWriter adiosWriter( ".", "channelFlow", storage );

   adiosWriter.add( pressure );
   adiosWriter.add( velocity );

   adiosWriter.write( maxLevel );
   /// [AdiosWriter]
#endif
}
