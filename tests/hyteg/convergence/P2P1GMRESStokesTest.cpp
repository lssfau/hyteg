#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "constant_stencil_operator/P2ConstantOperator.hpp"
#include "mixed_operator/ScalarToVectorOperator.hpp"
#include "mixed_operator/VectorLaplaceOperator.hpp"
#include "mixed_operator/VectorToScalarOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using namespace hyteg;

class CustomStokes : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ConstantVectorLaplaceOperator VelocityOperator_T;
   typedef P2ConstantVectorLaplaceOperator AOperatorType;
   typedef P2ToP1ConstantDivOperator       BOperatorType;
   typedef P1ToP2ConstantDivTOperator      BTOperatorType;
   typedef P1LumpedInvMassOperator         SchurOperatorType;

   CustomStokes( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , Lapl( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , invLumpedMass( storage, minLevel, maxLevel )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               const UpdateType                        updateType = Replace ) const
   {
      vertexdof::projectMean( src.p(), level );

      Lapl.apply( src.uvw(), dst.uvw(), level, flag, updateType );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );

      div.apply( src.uvw(), dst.p(), level, flag, updateType );

      vertexdof::projectMean( dst.p(), level );
   }

   const AOperatorType&     getA() const { return Lapl; }
   const BOperatorType&     getB() const { return div; }
   const BTOperatorType&    getBT() const { return divT; }
   const SchurOperatorType& getSchur() const { return invLumpedMass; }

 public:
   bool hasGlobalCells_;

   AOperatorType     Lapl;
   BOperatorType     div;
   BTOperatorType    divT;
   SchurOperatorType invLumpedMass;
};

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t minLevel_ = 0;
   const uint_t maxLevel_ = 0;
   const bool   blending  = false;

   // create the annulus mesh
   MeshInfo              meshInfo = MeshInfo::meshAnnulus( real_c( 1.0 ), real_c( 2.0 ), MeshInfo::meshFlavour::CRISS, 5, 3 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( blending )
   {
      AnnulusMap::setMap( setupStorage );
   }

   // set load balancing
   loadbalancing::roundRobinVolume( setupStorage );

   // set surface and cmb boundary flags by vertex location
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // create storage
   std::shared_ptr< hyteg::PrimitiveStorage > storage_ = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   // boundary conditions
   BoundaryCondition bc = BoundaryCondition::create0123BC();

   // fem functions
   P2P1TaylorHoodFunction< real_t > up_( "up", storage_, minLevel_, maxLevel_, bc );
   P2P1TaylorHoodFunction< real_t > f_( "f", storage_, minLevel_, maxLevel_, bc );
   P2P1TaylorHoodFunction< real_t > tmp_( "tmp", storage_, minLevel_, maxLevel_, bc );
   P2P1TaylorHoodFunction< real_t > refSol_( "refSol", storage_, minLevel_, maxLevel_, bc );

   // define functions
   std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return -2 * x[0] * x[1] - 2 * x[0] + std::pow( x[1], 2 ) + 23;
   };
   std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return 4 * std::pow( x[0], 2 ) - 8 * x[1] - 5;
   };
   std::function< real_t( const hyteg::Point3D& ) > uZ = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return 0;
   };

   std::function< real_t( const hyteg::Point3D& ) > p = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return x[0] - x[1] + 5;
   };

   std::function< real_t( const hyteg::Point3D& ) > ctrl_uX = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return -1;
   };

   std::function< real_t( const hyteg::Point3D& ) > ctrl_uY = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return -9;
   };

   std::function< real_t( const hyteg::Point3D& ) > ctrl_uZ = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return 0;
   };

   std::function< real_t( const hyteg::Point3D& ) > ctrl_p = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return 2 * x[1] + 10;
   };

   // operators
   auto stokesOperator_ = std::make_shared< CustomStokes >( storage_, minLevel_, maxLevel_ );
   auto P1Mass          = std::make_shared< P1ConstantMassOperator >( storage_, minLevel_, maxLevel_ );
   auto P2Mass          = std::make_shared< P2ConstantMassOperator >( storage_, minLevel_, maxLevel_ );
   auto Lapl_           = std::make_shared< CustomStokes::AOperatorType >( storage_, minLevel_, maxLevel_ );

   // interpolate
   {
      uint_t level = maxLevel_;

      tmp_.uvw()[0].interpolate( ctrl_uX, level, All );
      tmp_.uvw()[1].interpolate( ctrl_uY, level, All );
      if ( storage_->hasGlobalCells() )
      {
         tmp_.uvw()[2].interpolate( ctrl_uZ, level, All );
      }
      tmp_.p().interpolate( ctrl_p, level, All );

      up_.uvw()[0].interpolate( uX, level, DirichletBoundary );
      up_.uvw()[1].interpolate( uY, level, DirichletBoundary );
      if ( storage_->hasGlobalCells() )
      {
         up_.uvw()[2].interpolate( uZ, level, DirichletBoundary );
      }
      up_.p().interpolate( p, level, DirichletBoundary );

      refSol_.uvw()[0].interpolate( uX, level, All );
      refSol_.uvw()[1].interpolate( uY, level, All );
      if ( storage_->hasGlobalCells() )
      {
         refSol_.uvw()[2].interpolate( uZ, level, All );
      }
      refSol_.p().interpolate( p, level, All );

      P2Mass->apply( tmp_.uvw()[0], f_.uvw()[0], level, All );
      P2Mass->apply( tmp_.uvw()[1], f_.uvw()[1], level, All );
      if ( storage_->hasGlobalCells() )
      {
         P2Mass->apply( tmp_.uvw()[2], f_.uvw()[2], level, All );
      }
      P1Mass->apply( tmp_.p(), f_.p(), level, All );

      vertexdof::projectMean( refSol_.p(), level );
      vertexdof::projectMean( f_.p(), level );
   }

   tmp_.interpolate( 0, maxLevel_, All );

   // solver
   auto solver = hyteg::GMRESSolver< CustomStokes >( storage_, minLevel_, maxLevel_, 1000, 1000, 1e-16, 1e-16, 0 );
   //solver.setPrintInfo( true );

   // solve
   P2P1TaylorHoodFunction< real_t > Merr( "Merr", storage_, maxLevel_, maxLevel_ );
   DoFType                          flag = Inner | NeumannBoundary | FreeslipBoundary;

   solver.solve( *stokesOperator_, up_, f_, maxLevel_ );

   // check error
   tmp_.assign( { 1.0 }, { up_ }, maxLevel_, All );

   tmp_.assign( { 1.0, -1.0 }, { tmp_, refSol_ }, maxLevel_, All );

   for ( uint_t k = 0; k < tmp_.uvw().getDimension(); k++ )
   {
      P2Mass->apply( tmp_.uvw()[k], Merr.uvw()[k], maxLevel_, All );
   }
   P1Mass->apply( tmp_.p(), Merr.p(), maxLevel_, All );

   real_t discr_l2_err = std::sqrt( tmp_.dotGlobal( Merr, maxLevel_, flag ) );
   // real_t discr_l2_err_uvw = std::sqrt( tmp_.uvw().dotGlobal( Merr.uvw(), maxLevel_, flag ) );
   // real_t discr_l2_err_p   = std::sqrt( tmp_.p().dotGlobal( Merr.p(), maxLevel_, flag ) );

   // check residual
   stokesOperator_->apply( up_, tmp_, maxLevel_, flag );
   tmp_.assign( { 1.0, -1.0 }, { tmp_, f_ }, maxLevel_, flag );

   real_t stokesResidual  = std::sqrt( tmp_.dotGlobal( tmp_, maxLevel_, flag ) );
   real_t stokesResidualU = std::sqrt( tmp_.uvw().dotGlobal( tmp_.uvw(), maxLevel_, flag ) );
   real_t stokesResidualP = std::sqrt( tmp_.p().dotGlobal( tmp_.p(), maxLevel_, flag ) );

   // output
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "[Callback] abs. res: %10.5e | abs. res U: %10.5e | abs. res P: %10.5e | L2 Error: %10.5e",
                         stokesResidual,
                         stokesResidualU,
                         stokesResidualP,
                         discr_l2_err ) );

   // check error
   bool dp = std::is_same< real_t, double >();
   WALBERLA_CHECK_LESS( discr_l2_err, dp ? 1e-10 : 1e-3 );
   return EXIT_SUCCESS;
}
