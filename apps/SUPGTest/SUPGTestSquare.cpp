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
#include "hyteg_operators/operators/advection/P2ElementwiseAdvection.hpp"
#include "hyteg_operators/operators/supg_advection/P2ElementwiseSupgAdvection.hpp"
#include "hyteg_operators/operators/supg_diffusion/P2ElementwiseSupgDiffusion.hpp"

#include "mixed_operator/VectorMassOperator.hpp"

namespace hyteg {

class P2TransportOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
   using AdvectionOperatorType     = operatorgeneration::P2ElementwiseAdvection;
   using SUPGAdvectionOperatorType = operatorgeneration::P2ElementwiseSupgAdvection;
   using SUPGDiffusionOperatorType = operatorgeneration::P2ElementwiseSupgDiffusion;

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
          storage, minLevel, maxLevel, cpAdv, uVec.component( 0u ), uVec.component( 1u ) );

      cpAdvXdeltaP2.multElementwise( { cpAdv, deltaP2 }, maxLevel, All );
      kDiffXdeltaP2.multElementwise( { kDiff, deltaP2 }, maxLevel, All );

      advSUPG = std::make_shared< SUPGAdvectionOperatorType >(
          storage, minLevel, maxLevel, cpAdvXdeltaP2, uVec.component( 0u ), uVec.component( 1u ) );
      diffusionSUPG = std::make_shared< SUPGDiffusionOperatorType >(
          storage, minLevel, maxLevel, kDiffXdeltaP2, uVec.component( 0u ), uVec.component( 1u ) );
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
      cfg->readParameterFile( "./SUPGTestSquare.prm" );
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

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( -1.0, -1.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISSCROSS, 10u, 10u );
   // MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/SquareWithHole.msh" );
   // MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/SquareWithHoleBetterTriangles.msh" );
   // MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/SquareWithHoleFine.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // Point3D circleCenter{ { 0.0, 0.0, 0 } };
   // real_t  circleRadius = 0.3;

   // for ( const auto& it : setupStorage.getFaces() )
   // {
   //    Face& face = *( it.second );

   //    std::vector< PrimitiveID > neighborEdgesOnBoundary = face.neighborEdges();
   //    neighborEdgesOnBoundary.erase(
   //        std::remove_if( neighborEdgesOnBoundary.begin(),
   //                        neighborEdgesOnBoundary.end(),
   //                        [&setupStorage]( const PrimitiveID& id ) { return !setupStorage.onBoundary( id ); } ),
   //        neighborEdgesOnBoundary.end() );

   //    if ( neighborEdgesOnBoundary.size() > 0 )
   //    {
   //       Edge& edge = *setupStorage.getEdge( neighborEdgesOnBoundary[0] );

   //       const Vertex& vertex = *setupStorage.getVertex( face.get_vertex_opposite_to_edge( edge.getID() ) );

   //       if ( ( edge.getCoordinates()[0] - circleCenter ).norm() < 0.3 + 1e-4 )
   //       {
   //          // setupStorage.setGeometryMap(
   //          //     edge.getID(), std::make_shared< CircularMap >( edge.getCoordinates(), vertex.getCoordinates(), circleRadius ) );
   //          // setupStorage.setGeometryMap(
   //          //     face.getID(), std::make_shared< CircularMap >( edge.getCoordinates(), vertex.getCoordinates(), circleRadius ) );

   //          setupStorage.setGeometryMap( edge.getID(),
   //                                       std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
   //          setupStorage.setGeometryMap( face.getID(),
   //                                       std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
   //       }
   //    }
   // }

   const real_t tol = 1e-6;

   std::function< bool( const Point3D& ) > boundaryMarker = [tol]( const Point3D& x ) {
      if ( std::abs( x[0] + 1.0 ) < tol || std::abs( x[0] - 1.0 ) < tol || std::abs( x[1] + 1.0 ) < tol ||
           std::abs( x[1] - 1.0 ) < tol )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   setupStorage.setMeshBoundaryFlagsByCentroidLocation( DIRICHLETWALLS, boundaryMarker );

   // hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 3 );

   const uint_t level = mainConf->getParameter< uint_t >( "level" );
   const real_t eps   = 5e-3;

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
          real_t u = std::sqrt( v[0] * v[0] + v[1] * v[1] );

          return getDelta( u );
       };

   BoundaryCondition bcTemp;
   bcTemp.createDirichletBC( "DirichletWalls", DIRICHLETWALLS );

   P2Function< real_t > T( "T", storage, level, level, bcTemp );
   P2Function< real_t > fT( "fT", storage, level, level, bcTemp );

   P2VectorFunction< real_t > u( "u", storage, level, level );

   bool SUPG = false;

   if ( mainConf->getParameter< bool >( "SUPG" ) )
   {
      SUPG = true;
   }

   P2Function< real_t > cpAdv( "cpAdv", storage, level, level );
   P2Function< real_t > kDiff( "kDiff", storage, level, level );
   P2Function< real_t > deltaP2( "deltaP2", storage, level, level );

   cpAdv.interpolate( 1.0, level, All );
   kDiff.interpolate( eps, level, All );

   u.interpolate( { -std::sin( walberla::math::pi / 6 ), std::cos( walberla::math::pi / 6 ) }, level, All );

   deltaP2.interpolate( supgDelta, { u.component( 0u ), u.component( 1u ) }, level, All );

   P2TransportOperator p2Transport( storage, u, cpAdv, kDiff, deltaP2, level, level, eps, SUPG );

   std::function< real_t( const Point3D& ) > TDirichlet = []( const Point3D& x ) {
      if ( std::abs( x[0] - 1.0 ) < 1e-5 || ( std::abs( x[1] + 1.0 ) < 1e-5 && x[0] > 0.5 ) )
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

   GMRESSolver< P2TransportOperator > gmresTransport( storage, level, level, 2000UL, 2000UL, 1e-8, 1e-8 );

   gmresTransport.setPrintInfo( true );
   // p2Transport.apply(T, fT, level, All);
   gmresTransport.solve( p2Transport, T, fT, level );

   std::string vtkFilename = "SUPGTest";

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