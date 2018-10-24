
/**
 *
 * References:
 *
 * [1] Haines 2011: The Jeffery–Hamel similarity solution and its relation to flow in a diverging channel
 * [2] Fraenkel 1962: Laminar flow in symmetrical channels with slightly curved walls
 *                    I. On the Jeffery-Hamel solutions for flow between plane walls
 *
 */

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"

#include "tinyhhg_core/VTKWriter.hpp"

namespace hhg {

using walberla::real_c;
using walberla::real_t;
using walberla::math::PI;

void jefferyHamelFlowTest()
{
  ////////////////
  // Parameters //
  ////////////////

  // All variable names are defined as in [1].

  // free parameters

  const double eta      = 100.0;
  const double alpha    = 15.0 * (PI / 180.0);
  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  // derived parameters

  const double Rhat_o   = 1.0;
  const double Rhat_i   = Rhat_o / eta;
  const double epsilon  = 1e-8;

  ////////////
  // Domain //
  ////////////

  const auto meshInfo = MeshInfo::meshAnnulus( Rhat_i, Rhat_o, - alpha, alpha, MeshInfo::meshFlavour::CRISSCROSS, 4, 4 );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  // setting mesh boundary flags

  BoundaryCondition bc;
  const auto noSlipBCUID  = bc.createDirichletBC( "no-slip", 1 );
  const auto inflowBCUID  = bc.createDirichletBC( "inflow",  2 );
  const auto outflowBCUID = bc.createNeumannBC  ( "outflow", 3 );
  auto onInflowBoundary  = [ Rhat_i, epsilon ]( const Point3D & p ) { return std::abs( p.norm() - Rhat_i ) < epsilon; };
  auto onOutflowBoundary = [ Rhat_o, epsilon ]( const Point3D & p ) { return std::abs( p.norm() - Rhat_o ) < epsilon; };
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, onOutflowBoundary );
  setupStorage.setMeshBoundaryFlagsByVertexLocation( 3, onInflowBoundary );

  const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

  writeDomainPartitioningVTK( storage, "../../output", "JefferyHamel_Domain" );

  ///////////////
  // Functions //
  ///////////////

  P2P1TaylorHoodFunction< real_t > x            ( "x",             storage, minLevel, maxLevel );
  P2P1TaylorHoodFunction< real_t > exactSolution( "exactSolution", storage, minLevel, maxLevel );

  P2P1TaylorHoodStokesOperator stokesOperator( storage, minLevel, maxLevel );

  //////////////
  // Boundary //
  //////////////

  auto inflowProfile = [ alpha ]( const Point3D & p, const uint_t & coordinate ) -> real_t {
    const auto pSpherical = math::toSpherical( p );
    const real_t phiRescaled = p[2] / alpha;
    const real_t uSpherical = ( real_c( 3 ) / ( real_c( 4 ) * alpha ) ) * ( real_c( 1 ) - phiRescaled * phiRescaled );
    return math::toCartesian( Point3D( { uSpherical, real_c(0), real_c(0) } ) )[ coordinate ];
  };

  auto inflowProfileU = [ inflowProfile ]( const Point3D & p ) -> real_t { return inflowProfile( p, 0 ); };
  auto inflowProfileV = [ inflowProfile ]( const Point3D & p ) -> real_t { return inflowProfile( p, 1 ); };

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Interpolating inflow..." );

  x.u.interpolate( inflowProfileU, maxLevel, inflowBCUID );
  x.v.interpolate( inflowProfileV, maxLevel, inflowBCUID );

  /////////
  // VTK //
  /////////

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Setting up VTK..." );

  VTKOutput vtkOutput( "../../output", "JefferyHamel", storage );
  vtkOutput.add( &x.u );
  vtkOutput.add( &x.v );

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Writing initial VTK..." );

  vtkOutput.write( maxLevel, 0 );

}

}

int main()
{
  hhg::jefferyHamelFlowTest();
  return EXIT_SUCCESS;
}