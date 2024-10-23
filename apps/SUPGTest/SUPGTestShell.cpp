#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
// #include "hyteg/geometry/CircularMapNew.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"
#include "hyteg_operators/operators/advection/P2ElementwiseAdvectionIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/supg_advection/P2ElementwiseSupgAdvectionIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/supg_diffusion/P2ElementwiseSupgDiffusionIcosahedralShellMap.hpp"

#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/operators/P2TransportTALAOperator.hpp"

namespace hyteg {

class P2TransportOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
   using AdvectionOperatorType     = operatorgeneration::P2ElementwiseAdvectionIcosahedralShellMap;
   using SUPGAdvectionOperatorType = operatorgeneration::P2ElementwiseSupgAdvectionIcosahedralShellMap;
   using SUPGDiffusionOperatorType = operatorgeneration::P2ElementwiseSupgDiffusionIcosahedralShellMap;

 public:
   P2ElementwiseBlendingLaplaceOperator         lapl;
   std::shared_ptr< AdvectionOperatorType >     adv;
   std::shared_ptr< SUPGAdvectionOperatorType > advSUPG;
   std::shared_ptr< SUPGDiffusionOperatorType > diffusionSUPG;

   P2VectorFunction< real_t >& uVec;
   P2Function< real_t >&       cpAdv;
   P2Function< real_t >&       kDiff;
   P2Function< real_t >&       deltaP2;

   P2Function< real_t > cpAdvXdeltaP2;
   P2Function< real_t > kDiffXdeltaP2;

   const real_t k;
   const bool   supg;

   P2TransportOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                        P2VectorFunction< real_t >&                uVec_,
                        P2Function< real_t >&                      cpAdv_,
                        P2Function< real_t >&                      kDiff_,
                        P2Function< real_t >&                      deltaP2_,
                        size_t                                     minLevel,
                        size_t                                     maxLevel,
                        real_t                                     k_,
                        bool                                       supg_ )
   : Operator( storage, minLevel, maxLevel )
   , lapl( storage, minLevel, maxLevel )
   , uVec( uVec_ )
   , cpAdv( cpAdv_ )
   , kDiff( kDiff_ )
   , deltaP2( deltaP2_ )
   , cpAdvXdeltaP2( "cpAdvXdeltaP2", storage, minLevel, maxLevel )
   , kDiffXdeltaP2( "kDiffXdeltaP2", storage, minLevel, maxLevel )
   , k( k_ )
   , supg( supg_ )
   {
      adv = std::make_shared< AdvectionOperatorType >(
          storage, minLevel, maxLevel, cpAdv, uVec.component( 0u ), uVec.component( 1u ), uVec.component( 2u ) );

      cpAdvXdeltaP2.multElementwise( { cpAdv, deltaP2 }, maxLevel, All );
      kDiffXdeltaP2.multElementwise( { kDiff, deltaP2 }, maxLevel, All );

      advSUPG = std::make_shared< SUPGAdvectionOperatorType >(
          storage, minLevel, maxLevel, cpAdvXdeltaP2, uVec.component( 0u ), uVec.component( 1u ), uVec.component( 2u ) );
      diffusionSUPG = std::make_shared< SUPGDiffusionOperatorType >(
          storage, minLevel, maxLevel, kDiffXdeltaP2, uVec.component( 0u ), uVec.component( 1u ), uVec.component( 2u ) );
   }

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               const uint_t                level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      lapl.apply( src, dst, level, flag, updateType );

      dst.assign( { k }, { dst }, level, flag );

      if ( supg )
      {
         diffusionSUPG->apply( src, dst, level, flag, Add );
      }

      adv->apply( src, dst, level, flag, Add );

      if ( supg )
      {
         advSUPG->apply( src, dst, level, flag, Add );
      }
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  src,
                  const P2Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      auto matDiffusion = mat->createCopy();
      auto matAdvection = mat->createCopy();

      lapl.toMatrix( matDiffusion, src, dst, level, flag );

      if ( supg )
      {
         adv->toMatrix( matAdvection, src, dst, level, flag );
         advSUPG->toMatrix( matAdvection, src, dst, level, flag );
      }
      else
      {
         adv->toMatrix( matAdvection, src, dst, level, flag );
      }

      mat->createFromMatrixLinComb( { k, 1.0 }, { matDiffusion, matAdvection } );
   }
};

} // namespace hyteg

using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#if defined( HYTEG_BUILD_WITH_PETSC )
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./SUPGTestShell.prm" );
   }
   else
   {
      cfg = env.config();
   }

   enum BoundaryConditionsFlag
   {
      DIRICHLETWALLS = 17
   };

   auto mainConf = std::make_shared< walberla::Config::BlockHandle >( cfg->getBlock( "Parameters" ) );

   real_t rMin = 0.5, rMax = 1.0;

   MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 3u, 2u, rMin, rMax, MeshInfo::SHELLMESH_CLASSIC );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   IcosahedralShellMap::setMap( setupStorage );
   // hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 3 );

   const uint_t level = mainConf->getParameter< uint_t >( "level" );
   const real_t eps   = mainConf->getParameter< real_t >( "eps" );
   ;

   const real_t hk = MeshQuality::getMaximalEdgeLength( storage, level );

   std::function< real_t( real_t ) > getDelta = [hk, eps]( real_t v ) {
      real_t h = hk;
      real_t k = eps;

      real_t Pe = h * v / ( 2 * k );

      real_t xi;

      // replace xi with approximations in case Pe is too small or too large
      if ( Pe <= 0.5 )
      {
         // error smaller than ~1e-5 here
         xi = Pe / 3.0 - ( Pe * Pe * Pe ) / 45.0;
      }
      else if ( Pe >= 20.0 )
      {
         // error smaller than ~1e-15 here
         xi = 1.0 - 1.0 / Pe;
      }
      else
      {
         xi = 1.0 + 2.0 / ( std::exp( 2.0 * Pe ) - 1.0 ) - 1.0 / Pe;
      }

      real_t SUPG_scaling_ = 1.0;

      if ( v < 1e-6 )
      {
         return 0.0;
      }

      return SUPG_scaling_ * h / ( 2 * v ) * xi;
   };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > supgDelta =
       [getDelta]( const Point3D&, const std::vector< real_t >& v ) {
          real_t u = std::sqrt( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] );

          return getDelta( u );
       };

   BoundaryCondition bcTemp;
   bcTemp.createDirichletBC( "DirichletWalls", { MeshInfo::flagInnerBoundary, MeshInfo::flagOuterBoundary } );

   P2Function< real_t > T( "T", storage, level, level, bcTemp );
   P2Function< real_t > fT( "fT", storage, level, level, bcTemp );

   P2VectorFunction< real_t > u( "u", storage, level, level );

   bool SUPG = false;

   if ( mainConf->getParameter< bool >( "SUPG" ) )
   {
      SUPG = true;
   }

   auto cpAdv   = std::make_shared< P2Function< real_t > >( "cpAdv", storage, level, level );
   auto kDiff   = std::make_shared< P2Function< real_t > >( "kDiff", storage, level, level );
   auto deltaP2 = std::make_shared< P2Function< real_t > >( "deltaP2", storage, level, level );
   auto viscP2  = std::make_shared< P2Function< real_t > >( "viscP2", storage, level, level );

   auto uPtr = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPtr", storage, level, level );

   viscP2->interpolate( 1.0, level, All );
   cpAdv->interpolate( 1.0, level, All );
   kDiff->interpolate( eps, level, All );

   std::function< real_t( const Point3D& ) > uX = []( const Point3D& x ) { return x[0] / x.norm(); };

   std::function< real_t( const Point3D& ) > uY = []( const Point3D& x ) { return x[1] / x.norm(); };

   std::function< real_t( const Point3D& ) > uZ = []( const Point3D& x ) { return x[2] / x.norm(); };

   u.interpolate( { uX, uY, uZ }, level, All );

   uPtr->uvw().interpolate( { uX, uY, uZ }, level, All );

   deltaP2->interpolate( supgDelta, { u.component( 0u ), u.component( 1u ), u.component( 2u ) }, level, All );

   P2TransportOperator p2Transport( storage, u, *cpAdv, *kDiff, *deltaP2, level, level, eps, SUPG );

   terraneo::P2TransportIcosahedralShellMapOperator p2TransportFull( storage, level, level );

   p2TransportFull.setDiffusivityCoeff( kDiff );
   p2TransportFull.setAdvectionCoeff( cpAdv );
   p2TransportFull.setVelocity( uPtr );
   p2TransportFull.setViscosity( viscP2 );

   p2TransportFull.setTALADict( { { terraneo::TransportOperatorTermKey::SHEAR_HEATING_TERM, false },
                                  { terraneo::TransportOperatorTermKey::ADIABATIC_HEATING_TERM, false },
                                  { terraneo::TransportOperatorTermKey::INTERNAL_HEATING_TERM, false },
                                  { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_MMOC, false },
                                  { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, true },
                                  { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
                                  { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, SUPG } } );

   p2TransportFull.initializeOperators();

   p2TransportFull.setTimestep( 1.0 );

   std::function< real_t( const Point3D& ) > TDirichlet = [rMin, rMax]( const Point3D& x ) {
      real_t r = x.norm();
      if ( std::abs( r - rMax ) < 1e-5 )
      {
         return 0.0;
      }
      else if ( std::abs( r - rMin ) < 1e-5 )
      {
         return 1.0;
      }
      else
      {
         return 0.0;
      }
   };

   T.interpolate( TDirichlet, level, DirichletBoundary );

   fT.interpolate( 0.0, level, All );

   GMRESSolver< P2TransportOperator > gmresTransport( storage, level, level, 200UL, 200UL, 1e-8, 1e-8 );

   GMRESSolver< terraneo::P2TransportIcosahedralShellMapOperator > gmresFullTransport(
       storage, level, level, 200UL, 200UL, 1e-8, 1e-8 );

   gmresTransport.setPrintInfo( true );
   gmresFullTransport.setPrintInfo( true );
   // p2Transport.apply(T, fT, level, All);
   // gmresTransport.solve( p2Transport, T, fT, level );
   gmresFullTransport.solve( p2TransportFull, T, fT, level );

   std::string vtkFilename = "SUPGTestShell";

   if ( SUPG )
   {
      vtkFilename += "SUPG";
   }

   VTKOutput vtkOutput( "./output", vtkFilename, storage );

   vtkOutput.add( T );
   vtkOutput.add( fT );
   vtkOutput.add( u );

   vtkOutput.add( T.getVertexDoFFunction() );

   vtkOutput.write( level );

   return 0;
}