#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/solvers/Solver.hpp"
#include "tinyhhg_core/gridtransferoperators/ProlongationOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/RestrictionOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg
{

using walberla::uint_t;
using walberla::real_t;

template< class OperatorType >
class GeometricMultigridSolver
{
public:

  enum class CycleType
  {
    VCYCLE,
    WCYCLE
  };

  typedef typename OperatorType::srcType FunctionType;

//  static_assert( !std::is_same< FunctionType, typename OperatorType::dstType >::value,
//                 "CGSolver does not work for Operator with different src and dst FunctionTypes" );

  GeometricMultigridSolver( const std::shared_ptr< PrimitiveStorage >&              storage,
                            std::shared_ptr< Solver< OperatorType > >               smoother,
                            std::shared_ptr< Solver< OperatorType > >               coarseSolver,
                            std::shared_ptr< RestrictionOperator< FunctionType > >  restrictionOperator,
                            std::shared_ptr< ProlongationOperator< FunctionType > > prolongationOperator,
                            uint_t                                                  minLevel,
                            uint_t                                                  maxLevel,
                            uint_t                                                  nuPre  = 3,
                            uint_t                                                  nuPost = 3 )
  : minLevel_( minLevel )
  , maxLevel_( maxLevel )
  , smoother_( smoother )
  , coarseSolver_( coarseSolver )
  , restrictionOperator_( restrictionOperator )
  , prolongationOperator_( prolongationOperator )
  , ax_( "gmg_ax", storage, minLevel, maxLevel )
  , tmp_( "gmg_tmp", storage, minLevel, maxLevel )
  , nuPre_( nuPre )
  , nuPost_( nuPost )
  , flag_( hhg::Inner | hhg::NeumannBoundary)
  {
     zero_ = []( const hhg::Point3D& ) { return 0.0; };
  }

  ~GeometricMultigridSolver() = default;

  void solve(const OperatorType& A, FunctionType& x, FunctionType& b,const uint_t& level)
  {

    if (level == minLevel_)
    {
      coarseSolver_->solve(A, x, b, minLevel_);
    }
    else
    {
      // pre-smooth
      for (size_t i = 0; i < nuPre_; ++i)
      {
        smoother_->solve(A, x, b, level );
      }

      A.apply(x, ax_, level, flag_);
      tmp_.assign({1.0, -1.0}, { b, ax_ }, level, flag_);

      // restrict
      restrictionOperator_->restrict( tmp_, level, flag_ );

      b.assign({1.0}, { tmp_ }, level - 1, flag_);

      x.interpolate(zero_, level-1);

      solve(A, x, b, level-1);

      if (cycleType_ == CycleType::WCYCLE) {
        solve(A, x, b, level-1);
      }

      // prolongate
      tmp_.assign({1.0}, { x }, level, flag_);
      prolongationOperator_->prolongate( x, level-1, flag_ );
      x.add({1.0}, { tmp_ }, level, flag_);

      // post-smooth
      for (size_t i = 0; i < nuPost_; ++i)
      {
        smoother_->solve(A, x, b, level );
      }
    }

  }

private:

  uint_t nuPre_;
  uint_t nuPost_;
  uint_t minLevel_;
  uint_t maxLevel_;

  hhg::DoFType flag_;
  CycleType cycleType_;

  std::shared_ptr< hhg::Solver< OperatorType > > smoother_;
  std::shared_ptr< hhg::Solver< OperatorType > > coarseSolver_;
  std::shared_ptr< hhg::RestrictionOperator< FunctionType > > restrictionOperator_;
  std::shared_ptr< hhg::ProlongationOperator< FunctionType > >  prolongationOperator_;

  FunctionType ax_;
  FunctionType tmp_;

  std::function<real_t(const hhg::Point3D&)> zero_;

};

}
