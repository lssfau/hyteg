
#pragma once

#include "tinyhhg_core/solvers/Solver.hpp"

namespace hyteg {

template < typename OperatorType >
class SolverLoop : public Solver< OperatorType >
{
public:

    typedef typename OperatorType::srcType FunctionType;

    SolverLoop( const std::shared_ptr< Solver< OperatorType > > & solver,
                const uint_t & iterations ) :
                solver_( solver ),
                iterations_( iterations )
    {}

    void solve( const OperatorType& A, const FunctionType& x, const FunctionType& b, const uint_t level ) override
    {
      for ( uint_t i = 0; i < iterations_; i++ )
        solver_->solve( A, x, b, level );
    }

private:

    std::shared_ptr< Solver< OperatorType > > solver_;
    uint_t iterations_;
};

}