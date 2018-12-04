#pragma once

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/StokesOperatorTraits.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/ProlongationOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/RestrictionOperator.hpp"

#include "MinresSolver.hpp"

namespace hhg {

template < class Function, class Operator, class CoarseGridSolver, bool Tensor >
class UzawaSolver
{
 public:
   enum class CycleType
   {
      VCYCLE,
      WCYCLE
   };

   UzawaSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                const CoarseGridSolver&                    coarseGridSolver,
                const uint_t&                              minLevel,
                const uint_t&                              maxLevel,
                const uint_t&                              numberOfPreSmoothingSteps,
                const uint_t&                              numberOfPostSmoothingSteps,
                const uint_t&                              smoothingStepIncrement,
                const real_t&                              relaxParam = walberla::real_c( 0.3 ) )
   : coarseGridSolver_( coarseGridSolver )
   , pspg_( storage, minLevel, maxLevel )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , nuPre_( numberOfPreSmoothingSteps )
   , nuPost_( numberOfPostSmoothingSteps )
   , nuAdd_( smoothingStepIncrement )
   , ax_( "uzw_ax", storage, minLevel, maxLevel )
   , tmp_( "uzw_tmp", storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , relaxParam_( relaxParam )
   {
      zero_ = []( const hhg::Point3D& ) { return 0.0; };
      if( std::is_same< Function, P2P1TaylorHoodFunction< real_t > >::value )
      {
         restrictionOperator_  = std::make_shared< hhg::P2P1StokesToP2P1StokesRestriction >();
         prolongationOperator_ = std::make_shared< hhg::P2P1StokesToP2P1StokesProlongation >();
         prolongateFunction    = hhg::prolongateP2P1;
      } else if( std::is_same< Function, P1StokesFunction< real_t > >::value )
      {
         //prolongateFunction = hhg::prolongateP1;

         //        restrictionOperator_  = std::make_shared< hhg::P1P1StokesToP1P1StokesRestriction >();
         //        prolongationOperator_ = std::make_shared< hhg::P1P1StokesToP1P1StokesProlongation >();
      }
   }

   ~UzawaSolver() = default;

   void init( Function ) { WALBERLA_ABORT( " provide init function for discretization " ); };

   void solve( Operator& A,
               Function& x,
               Function& b,
               Function& r,
               uint_t    level,
               real_t    tolerance,
               size_t    maxiter,
               DoFType   flag      = All,
               CycleType cycleType = CycleType::VCYCLE,
               bool      printInfo = false )
   {
      if( level == minLevel_ )
      {
         coarseGridSolver_.solve( A, x, b, r, level, 1e-16, maxiter, flag, printInfo );
         //      uzawaSmooth(A, x, b, r, level, flag);
      } else
      {
         // pre-smooth
         for( size_t i = 0; i < nuPre_; ++i )
         {
            uzawaSmooth( A, x, b, r, level, flag );
         }

         A.apply( x, ax_, level, flag );
         r.assign( {1.0, -1.0}, {&b, &ax_}, level, flag );

         // restrict
         restrictionOperator_->restrict( r, level, flag );

         b.assign( {1.0}, {&r}, level - 1, flag );
         //      vertexdof::projectMean(b.p, ax_.p, level-1);

         x.interpolate( zero_, level - 1 );

         // solve on coarser level
         nuPre_ += nuAdd_;
         nuPost_ += nuAdd_;
         solve( A, x, b, r, level - 1, tolerance, maxiter, flag, cycleType, printInfo );
         if( cycleType == CycleType::WCYCLE )
         {
            solve( A, x, b, r, level - 1, tolerance, maxiter, flag, cycleType, printInfo );
         }
         nuPre_ -= nuAdd_;
         nuPost_ -= nuAdd_;

         // prolongate
         tmp_.assign( {1.0}, {&x}, level, flag );
         prolongationOperator_->prolongate( x, level - 1, flag );
         x.add( {1.0}, {&tmp_}, level, flag );

         // post-smooth
         for( size_t i = 0; i < nuPost_; ++i )
         {
            uzawaSmooth( A, x, b, r, level, flag );
         }
      }
   }

 private:
   void uzawaSmooth( Operator& A, Function& x, Function& b, Function& r, uint_t level, DoFType flag )
   {
      uzawaSmooth( A,
                   x,
                   b,
                   r,
                   level,
                   flag,
                   std::integral_constant< bool, Tensor >(),
                   std::integral_constant< bool, has_pspg_block< Operator >::value >() );
   }

   // Block-Laplace variant
   void uzawaSmooth( Operator& A,
                     Function& x,
                     Function& b,
                     Function& r,
                     uint_t    level,
                     DoFType   flag,
                     std::false_type /* tensor */,
                     std::true_type /* PSPG */ )
   {
      A.divT_x.apply( x.p, r.u, level, flag, Replace );
      r.u.assign( {1.0, -1.0}, {&b.u, &r.u}, level, flag );
      A.A.smooth_gs( x.u, r.u, level, flag );

      A.divT_y.apply( x.p, r.v, level, flag, Replace );
      r.v.assign( {1.0, -1.0}, {&b.v, &r.v}, level, flag );
      A.A.smooth_gs( x.v, r.v, level, flag );

      if( hasGlobalCells_ )
      {
         A.divT_z.apply( x.p, r.w, level, flag, Replace );
         r.w.assign( {1.0, -1.0}, {&b.w, &r.w}, level, flag );
         A.A.smooth_gs( x.w, r.w, level, flag );
      }

      A.div_x.apply( x.u, r.p, level, flag, Replace );
      A.div_y.apply( x.v, r.p, level, flag, Add );

      if( hasGlobalCells_ )
      {
         A.div_z.apply( x.w, r.p, level, flag, Add );
      }

      r.p.assign( {1.0, -1.0}, {&b.p, &r.p}, level, flag );

      A.pspg.smooth_sor( x.p, r.p, relaxParam_, level, flag );
   }

   // Tensor variant
   void uzawaSmooth( Operator& A,
                     Function& x,
                     Function& b,
                     Function& r,
                     uint_t    level,
                     DoFType   flag,
                     std::true_type /* tensor */,
                     std::true_type /* PSPG */ )
   {
      A.divT_x.apply( x.p, r.u, level, flag, Replace );
      A.A_uv.apply( x.v, r.u, level, flag, Add );
      r.u.assign( {1.0, -1.0}, {&b.u, &r.u}, level, flag );
      A.A_uu.smooth_gs( x.u, r.u, level, flag );

      A.divT_y.apply( x.p, r.v, level, flag, Replace );
      A.A_vu.apply( x.u, r.v, level, flag, Add );
      r.v.assign( {1.0, -1.0}, {&b.v, &r.v}, level, flag );
      A.A_vv.smooth_gs( x.v, r.v, level, flag );

      A.div_x.apply( x.u, r.p, level, flag, Replace );
      A.div_y.apply( x.v, r.p, level, flag, Add );

      r.p.assign( {1.0, -1.0}, {&b.p, &r.p}, level, flag );

      A.pspg.smooth_sor( x.p, r.p, relaxParam_, level, flag );
   }

   // Block-Laplace variant without stabilization
   void uzawaSmooth( Operator& A,
                     Function& x,
                     Function& b,
                     Function& r,
                     uint_t    level,
                     DoFType   flag,
                     std::false_type /* tensor */,
                     std::false_type /* PSPG */ )
   {
      A.divT_x.apply( x.p, r.u, level, flag, Replace );
      r.u.assign( {1.0, -1.0}, {&b.u, &r.u}, level, flag );
      A.A.smooth_gs( x.u, r.u, level, flag );

      A.divT_y.apply( x.p, r.v, level, flag, Replace );
      r.v.assign( {1.0, -1.0}, {&b.v, &r.v}, level, flag );
      A.A.smooth_gs( x.v, r.v, level, flag );

      A.divT_x.apply( x.p, r.u, level, flag, Replace );
      r.u.assign( {1.0, -1.0}, {&b.u, &r.u}, level, flag );
      A.A.smooth_gs( x.u, r.u, level, flag );

      A.divT_y.apply( x.p, r.v, level, flag, Replace );
      r.v.assign( {1.0, -1.0}, {&b.v, &r.v}, level, flag );
      A.A.smooth_gs( x.v, r.v, level, flag );

      A.div_x.apply( x.u, r.p, level, flag | DirichletBoundary, Replace );
      A.div_y.apply( x.v, r.p, level, flag | DirichletBoundary, Add );
      r.p.assign( {1.0, -1.0}, {&b.p, &r.p}, level, flag | DirichletBoundary );

      //    invPressureMass_.apply(r.p, tmp_.p, level, flag | DirichletBoundary);
      //    x.p.add({-relaxParam_}, {&tmp_.p}, level, flag | DirichletBoundary);

      //    tmp_.p.interpolate(0.0, level);
      //    pressureMass_.smooth_sor(tmp_.p, r.p, relaxParam_, level, flag | DirichletBoundary);
      //    x.p.add({-1.0}, {&tmp_.p}, level, flag | DirichletBoundary);

      tmp_.p.interpolate( 0.0, level );
      pspg_.smooth_sor( tmp_.p, r.p, relaxParam_, level, flag | DirichletBoundary );
      x.p.add( {1.0}, {&tmp_.p}, level, flag | DirichletBoundary );
   }

   uint_t nuPre_;
   uint_t nuPost_;
   uint_t nuAdd_;

   uint_t minLevel_;
   uint_t maxLevel_;

   real_t relaxParam_;

   CoarseGridSolver                                         coarseGridSolver_;
   std::shared_ptr< hhg::RestrictionOperator< Function > >  restrictionOperator_;
   std::shared_ptr< hhg::ProlongationOperator< Function > > prolongationOperator_;

   P1PSPGOperator pspg_;

   Function ax_;
   Function tmp_;

   bool hasGlobalCells_;

   std::function< real_t( const hhg::Point3D& ) > zero_;
};

//template<class Operator, class CoarseGridSolver, bool Tensor>
//class UzawaSolver<P2P1TaylorHoodFunction< real_t >, Operator, CoarseGridSolver, Tensor>::init(){
//    restrictionOperator_  = std::make_shared< hhg::P2P1StokesToP2P1StokesRestriction >();
//    prolongationOperator_ = std::make_shared< hhg::P2P1StokesToP2P1StokesProlongation >();
//  }

} // namespace hhg
