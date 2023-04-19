#include <complex>
#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"

#include "hyteg/composites/P1P0StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/egfunctionspace/EGConvTestUtils.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVersion.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

/* Multi sinker benchmark described in section 2.1 of [1]. Specify solver type, 
   number of elements in each spatial direction, number of sinkers, viscosity 
   contrast, decay rate of the sinkers and radius of the sinkers.

   [1]: "WEIGHTED BFBT PRECONDITIONER FOR STOKES FLOW PROBLEMS WITH HIGHLY HETEROGENEOUS VISCOSITY" 
   by JOHANN RUDI, GEORG STADLER, AND OMAR GHATTAS 
*/

template < typename StokesOperatorType >
void MultiSinker( const std::string& name,
                  const uint_t&      level,
                  const uint_t&      nxy,
                  const uint_t&      nSinkers,
                  const real_t&      DR,
                  const real_t&      delta,
                  const real_t&      omega )
{
   using StokesFunctionType          = typename StokesOperatorType::srcType;
   using StokesFunctionNumeratorType = typename StokesFunctionType::template FunctionType< idx_t >;
   real_t visc_min                   = std::pow( DR, -0.5 );
   real_t visc_max                   = std::pow( DR, 0.5 );

   // storage and domain
   auto meshInfo = MeshInfo::meshCuboid(Point3D({0, 0, 0}), Point3D({1, 1, 1}), nxy, nxy, nxy);
   //auto meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // some general info about the testcase
   StokesFunctionNumeratorType numerator( "numerator", storage, level, level + 1 );
   numerator.enumerate( level ); numerator.enumerate( level + 1);
   uint_t globalDoFs      = numberOfGlobalDoFs( numerator, level );
   uint_t globalDoFsFiner = numberOfGlobalDoFs( numerator, level + 1 );
   WALBERLA_LOG_INFO_ON_ROOT( "### Computing " << name << " on level " << level << " ###" );
   WALBERLA_LOG_INFO_ON_ROOT( "Global DoFs: " << globalDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "Global DoFs finer: " << globalDoFsFiner );
   WALBERLA_LOG_INFO_ON_ROOT( "#Sinkers: " << nSinkers << ", visc_max: " << visc_max );

   // function setup
   StokesFunctionType x( "x", storage, level, level + 1 );
   StokesFunctionType xFiner( "xFiner", storage, level, level + 1 );

   StokesFunctionType btmp( "btmp", storage, level, level + 1 );
   StokesFunctionType b( "b", storage, level, level + 1 );
   StokesFunctionType residuum( "res", storage, level, level + 1 );
   StokesFunctionType errEq( "errEq", storage, level, level + 1 );

   P2Function< real_t > viscFunc( "viscFunc", storage, level, level + 1 );
   using hyteg::dg::eg::copyBdry;
    if constexpr (hyteg::dg::eg::isEGP0Discr<StokesOperatorType>()) {
        copyBdry(x);
        copyBdry(xFiner);
        copyBdry(btmp);
        copyBdry(b);
        copyBdry(residuum);
        copyBdry(errEq);
    }

   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return real_c( 0 ); };
   x.uvw().interpolate( zero, level );
   x.p().interpolate( zero, level );

   // generate sinker centers randomly (without them exiting the domain)
   std::uniform_real_distribution< real_t > unif( 0.0, 1.0 );
   std::default_random_engine               re;
   re.seed( 25151 );

   std::vector< Point3D > centers;
   for ( uint_t c = 0; c < nSinkers; c++ )
   {
      centers.push_back( Point3D( { unif( re ), unif( re ), unif( re ) } ) );
   }

   // characteristic function for sinkers
   std::function< real_t( const hyteg::Point3D& ) > Xi = [centers, delta, omega]( const hyteg::Point3D& xx ) {
      real_t val = 1;
      for ( auto& c : centers )
      {
         auto distance = c - xx;
         val *= ( 1 - exp( -delta * std::pow( std::max( 0.0, distance.norm() - omega / 2 ), 2 ) ) );
      }
      return val;
   };

   // viscosity function
   std::function< real_t( const hyteg::Point3D& ) > viscosity = [visc_max, visc_min, Xi]( const hyteg::Point3D& xx ) {
      return ( visc_max - visc_min ) * ( 1 - Xi( xx ) ) + visc_min;
   };

   // right hand side: "force sinkers downward": negative v velocity at sinker locations
   std::function< real_t( const hyteg::Point3D& ) > rhsV = [Xi]( const hyteg::Point3D& xx ) { return 10 * ( Xi( xx ) - 1 ); };
   btmp.uvw().interpolate( { zero, zero, rhsV }, level, All );
   btmp.uvw().interpolate( { zero, zero, rhsV }, level + 1, All );

   if constexpr ( hyteg::dg::eg::isP2P1Discr< StokesOperatorType >() )
   {
      P2ConstantMassOperator M_vel( storage, level, level + 1 );
      M_vel.apply( btmp.uvw()[2], b.uvw()[2], level, All );
      M_vel.apply( btmp.uvw()[2], b.uvw()[2], level + 1, All );
   }
   else if constexpr ( hyteg::dg::eg::isEGP0Discr< StokesOperatorType >() )
   {
      hyteg::dg::eg::EGMassOperator M_vel( storage, level, level + 1 );
      M_vel.apply( btmp.uvw(), b.uvw(), level, All, Replace );
      M_vel.apply( btmp.uvw(), b.uvw(), level + 1, All, Replace );
   }

   // operator setup
   StokesOperatorType Op( storage, level, level, viscosity );

   StokesOperatorType OpFiner( storage, level + 1, level + 1, viscosity );
   // Visualization
   /* VTKOutput vtkOutput( "../../output", name, storage );
   viscFunc.interpolate( viscosity, level, All );
   vtkOutput.add( viscFunc );
   vtkOutput.add( x.uvw() );
   vtkOutput.add( x.p() );
   vtkOutput.add( xFiner.uvw() );
   vtkOutput.add( xFiner.p() );
   vtkOutput.add( b.uvw() );
   vtkOutput.add( b.p() );
   if constexpr ( hyteg::dg::eg::isEGP0Discr< StokesOperatorType >() )
   {
      vtkOutput.add( *x.uvw().getConformingPart() );
      vtkOutput.add( *x.uvw().getDiscontinuousPart() );
      vtkOutput.add( *xFiner.uvw().getConformingPart() );
      vtkOutput.add( *xFiner.uvw().getDiscontinuousPart() );
   }
*/

   if ( true )
   {
      {
         x.interpolate([&unif, &re](const hyteg::Point3D &) { return unif(re); }, level, hyteg::Inner);
         PETScBlockPreconditionedStokesSolver< StokesOperatorType > solver(
             storage, level, 1e-3, std::numeric_limits< PetscInt >::max(), 6, 1 );
         solver.disableApplicationBC( std::is_same< StokesOperatorType, dg::eg::EGP0EpsilonOperatorStokesNitscheBC >::value );
         solver.solve( Op, x, b, level );
      }
      WALBERLA_LOG_INFO_ON_ROOT( "First solve done." );

      {
         xFiner.interpolate([&unif, &re](const hyteg::Point3D &) { return unif(re); }, level, hyteg::Inner);
         PETScBlockPreconditionedStokesSolver< StokesOperatorType > solverFiner(
             storage, level + 1, 1e-5, std::numeric_limits< PetscInt >::max(), 6, 1 );
         solverFiner.disableApplicationBC(
             std::is_same< StokesOperatorType, dg::eg::EGP0EpsilonOperatorStokesNitscheBC >::value );
         solverFiner.solve( OpFiner, xFiner, b, level + 1 );
      }
   }

   // subtract mean
   if constexpr ( hyteg::dg::eg::isEGP0Discr< StokesOperatorType >() )
   {
      hyteg::dg::projectMean( x.p(), level );
      hyteg::dg::projectMean( xFiner.p(), level + 1 );
   }
   else if constexpr ( hyteg::dg::eg::isP2P1Discr< StokesOperatorType >() )
   {
      hyteg::vertexdof::projectMean( x.p(), level );
      hyteg::vertexdof::projectMean( xFiner.p(), level + 1 );
   }

   // error equivalent
   if constexpr ( hyteg::dg::eg::isP2P1Discr< StokesOperatorType >() )
   {
      P2toP2QuadraticProlongation P2P2ProlongationOp;
      for ( uint_t k = 0; k < x.uvw().getDimension(); k++ )
      {
         P2P2ProlongationOp.prolongate( x.uvw()[k], level, All );
      }
      errEq.assign( { 1.0, -1.0 }, { xFiner, x }, level + 1, Inner );

      WALBERLA_LOG_INFO_ON_ROOT( "error equivalent u = " << errEq.uvw().getMaxComponentMagnitude( level + 1, All ) );
   }
   else if constexpr ( hyteg::dg::eg::isEGP0Discr< StokesOperatorType >() )
   {
      P1toP1LinearProlongation P1P1ProlongationOp;
      for ( uint_t k = 0; k < x.uvw().getConformingPart()->getDimension(); k++ )
      {
         P1P1ProlongationOp.prolongate( ( *x.uvw().getConformingPart() )[k], level, All );
      }
      errEq.uvw().getConformingPart()->assign(
          { 1.0, -1.0 }, { *xFiner.uvw().getConformingPart(), *x.uvw().getConformingPart() }, level + 1, All );

      WALBERLA_LOG_INFO_ON_ROOT( "error equivalent u = " << errEq.uvw().getMaxMagnitude( level + 1 ) );
   }

   //vtkOutput.write( level, 1 );
}

} // namespace hyteg

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   uint_t nSinkers = 1;
   if ( true )
   {
      MultiSinker< hyteg::dg::eg::EGP0EpsilonOperatorStokesNitscheBC >( "MultiSinker_EGP0", 3, 1, nSinkers, 1000, 200, 0.1 );
   }

   if ( false )
   {
      MultiSinker< P2P1ElementwiseAffineEpsilonStokesOperator >( "MultiSinker_P2P1", 3, 1, nSinkers, 1000, 200, 0.1 );
   }

   return EXIT_SUCCESS;
}
