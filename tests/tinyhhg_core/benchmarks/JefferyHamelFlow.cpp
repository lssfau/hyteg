
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

  const double alpha    = 15.0 * (PI / 180.0);
  const double Rhat_i   = 0.3;
  const double Rhat_o   = 1.0;
  const double eta      = Rhat_o / Rhat_i;
  const double epsilon  = 1e-8;

  const uint_t minLevel = 2;
  const uint_t maxLevel = 4;

  ////////////
  // Domain //
  ////////////

  const auto meshInfo = MeshInfo::meshAnnulus( Rhat_i, Rhat_o, - alpha, alpha, MeshInfo::meshFlavour::CRISSCROSS, 4, 4 );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  // setting mesh boundary flags
  const auto bc = BoundaryCondition::create012BC();
  auto onOutflowBoundary = [ Rhat_o, epsilon ]( const Point3D & p ) { return std::abs( p.norm() - Rhat_o ) < epsilon; };
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, onOutflowBoundary );

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



}

}

int main()
{
  hhg::jefferyHamelFlowTest();
  return EXIT_SUCCESS;
}