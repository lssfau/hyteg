

#pragma once
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_invk_mass_affine_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_sqrtk_mass_affine_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_sqrtk_mass_affine_q6.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"

#include "mixed_operator/P1P1StokesOperator.hpp"
#include "mixed_operator/P2P1TaylorHoodStokesOperator.hpp"
//#include "hyteg/operators/combinedoperators/BFBTOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
// PETSc
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"

namespace hyteg {

/// \brief Contains convenience functions that assemble solvers for common problems.
namespace solvertemplates {

/// \brief Returns a pressure preconditioned MINRES solver for the Stokes system.
///
/// The pressure is pre-multiplied with the inverse of the lumped mass matrix.
/// It is assumed that the pressure is discretized with P1 finite elements.
///
/// \tparam StokesOperatorType most types of Stokes operators should be possible to pass here
///                            since the MINRES solver does not require any knowledge about the
///                            structure of the A-block
/// \param storage the PrimitiveStorage that defines the domain
/// \param level the refinement level of the grid
/// \param absoluteTargetResidual absolute (as opposed to relative) residual as a stopping criterion for the iteration
/// \param maxIterations if not converged to the target residual, the iteration stops after this many iterations
///
template < typename StokesOperatorType >
std::shared_ptr< Solver< StokesOperatorType > > stokesMinResSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                    const uint_t&                              level,
                                                                    const real_t& absoluteTargetResidual,
                                                                    const uint_t& maxIterations,
                                                                    bool          printInfo = false )
{
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< StokesOperatorType, P1LumpedInvMassOperator > >(
           storage, level, level );
   auto pressurePreconditionedMinResSolver = std::make_shared< MinResSolver< StokesOperatorType > >(
       storage, level, level, maxIterations, absoluteTargetResidual, pressurePreconditioner );

   pressurePreconditionedMinResSolver->setPrintInfo( printInfo );

   return pressurePreconditionedMinResSolver;
}

/// \brief Returns a block preconditioned MINRES (viscosity weighted mass matrix preconditioning for p) solver for the Stokes system with varying viscosity
///
///
/// The pressure is pre-multiplied with the inverse of the lumped mass matrix, weighted by the viscosity.
/// It is assumed that the pressure is discretized with P1 finite elements.
/// The Epsilon operator is inverted approximately via GMG cycle(s).
///
///
/// \param StokesOperatorType assumed to be
/// \param storage the PrimitiveStorage that defines the domain
/// \param level the refinement level of the grid
/// \param absoluteTargetResidual absolute (as opposed to relative) residual as a stopping criterion for the iteration
/// \param maxIterations if not converged to the target residual, the iteration stops after this many iterations
///
/*
template < typename StokesOperatorType >
std::shared_ptr< Solver< StokesOperatorType > >
    varViscStokesMinResSolver( const std::shared_ptr< PrimitiveStorage >&       storage,
                               const uint_t&                                    maxLevel,
                               std::function< real_t( const hyteg::Point3D& ) > viscosity,
                               const uint_t&                                    nPrecCycles,
                               const real_t&                                    relativeResidual,
                               const uint_t&                                    maxIterations,
                               bool                                             printInfo )
{
   using StokesFunctionType          = typename StokesOperatorType::srcType;
   using StokesFunctionNumeratorType = typename StokesFunctionType::template FunctionType< idx_t >;
   StokesFunctionNumeratorType Numerator( "Num", storage, maxLevel, maxLevel );
   Numerator.enumerate( maxLevel );

   auto LU = std::make_shared< PETScLUSolver< StokesOperatorType > >( storage, maxLevel );

   // construct pressure preconditioning operator: inverse, lumped, viscosity weighted, P1 mass matrix
   auto pPrecOp = std::make_shared< P1BlendingLumpedInverseDiagonalOperator >(
       storage,
       maxLevel,
       maxLevel,
       std::make_shared< P1RowSumForm >( std::make_shared< forms::p1_invk_mass_affine_q4 >( viscosity, viscosity ) ) );

   auto prec =
       std::make_shared< StokesBlockDiagonalPreconditioner< StokesOperatorType, P1BlendingLumpedInverseDiagonalOperator > >(
           storage, maxLevel, maxLevel, nPrecCycles, pPrecOp, LU );

   auto solver = std::make_shared< MinResSolver< StokesOperatorType > >(
       storage, maxLevel, maxLevel, maxIterations, relativeResidual, prec );
   solver->setPrintInfo( printInfo );
   return solver;
}
*/
/// \brief Returns a block preconditioned MINRES (weighted BFBT preconditioning for p) solver for the Stokes system with varying viscosity.
///
/// The pressure is pre-multiplied with the inverse of the wBFBT operator, a better approximation to the Schur complement than the viscosity weighted mass matrix.
/// It is assumed that the pressure is discretized with P1 finite elements.
/// The Epsilon operator is inverted approximately via GMG cycle(s).
///
///
/// \param StokesOperatorType assumed to be
/// \param storage the PrimitiveStorage that defines the domain
/// \param level the refinement level of the grid
/// \param absoluteTargetResidual absolute (as opposed to relative) residual as a stopping criterion for the iteration
/// \param maxIterations if not converged to the target residual, the iteration stops after this many iterations
///
/*
template < typename StokesOperatorType >
std::shared_ptr< Solver< StokesOperatorType > >
    BFBTStokesMinResSolver( const std::shared_ptr< PrimitiveStorage >&       storage,
                            const uint_t&                                    maxLevel,
                            std::function< real_t( const hyteg::Point3D& ) > viscosity,
                            const real_t&                                    relativeTolerance,
                            const uint_t&                                    maxIterations,
                            bool                                             printInfo,
                            const std::vector< BoundaryCondition >&          VelocitySpaceBCs

    )
{
   using StokesFunctionType          = typename StokesOperatorType::srcType;
   using StokesFunctionNumeratorType = typename StokesFunctionType::template FunctionType< idx_t >;
   StokesFunctionNumeratorType Numerator( "Num", storage, maxLevel, maxLevel );
   Numerator.enumerate( maxLevel );

   auto LU = std::make_shared< PETScLUSolver< StokesOperatorType > >( storage, maxLevel, Numerator );

   // construct pressure preconditioning operator: inverse, lumped, viscosity weighted, P2 mass matrix
   auto pPrecOp = std::make_shared< BFBT_P2P1 >( storage, maxLevel, maxLevel, viscosity, VelocitySpaceBCs );
   auto prec    = std::make_shared< StokesBlockDiagonalPreconditioner< StokesOperatorType, BFBT_P2P1 > >(
       storage, maxLevel, maxLevel, 1, pPrecOp, LU );

   auto solver = std::make_shared< MinResSolver< StokesOperatorType > >(
       storage, maxLevel, maxLevel, maxIterations, relativeTolerance, prec );
   solver->setPrintInfo( printInfo );
   return solver;
}
*/
/// \brief Returns a block preconditioned MINRES solver for the Stokes system with varying viscosity.
///
/// The pressure is pre-multiplied with the inverse of the lumped mass matrix.
/// It is assumed that the pressure is discretized with P1 finite elements.
/// The Epsilon operator is inverted approximately via GMG cycle(s).
///
///
/// \param StokesOperatorType assumed to be
/// \param storage the PrimitiveStorage that defines the domain
/// \param level the refinement level of the grid
/// \param absoluteTargetResidual absolute (as opposed to relative) residual as a stopping criterion for the iteration
/// \param maxIterations if not converged to the target residual, the iteration stops after this many iterations
///

/*
template<typename StokesOperatorType>
std::shared_ptr< Solver< StokesOperatorType > >
    blkdiagPrecStokesMinResSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                   const uint_t&                              minLevel,
                                   const uint_t&                              maxLevel,
                                   const real_t&                              absoluteTargetResidual,
                                   const real_t&                              relativeVelBlockResidual,
                                   const uint_t&                              maxIterations,
                                   //const real_t&                                        relax,
                                   bool printInfo )
{
   // Velocity block solver
   //auto CG = std::make_shared<PETScCGSolver< P2ElementwiseAffineEpsilonOperator > >( storage, maxLevel, relativeVelBlockResidual, 1e-12);
   auto LU               = std::make_shared< PETScLUSolver< P2ElementwiseAffineEpsilonOperator > >( storage, maxLevel );
   auto coarseGridSolver = std::make_shared< PETScLUSolver< hyteg::P2ElementwiseAffineEpsilonOperator > >( storage, minLevel );
   auto smoother         = std::make_shared< WeightedJacobiSmoother< hyteg::P2ElementwiseAffineEpsilonOperator > >(
       storage, minLevel, maxLevel, 1.0 );
   auto prolongationOperator = std::make_shared< P2toP2VectorProlongation >();
   auto restrictionOperator  = std::make_shared< P2toP2VectorRestriction >();
   auto gmgSolver            = std::make_shared< GeometricMultigridSolver< hyteg::P2ElementwiseAffineEpsilonOperator > >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3, 3 );
   auto CG = std::make_shared< CGSolver< P2ElementwiseAffineEpsilonOperator > >(
       storage, minLevel, maxLevel, std::numeric_limits< uint_t >::max(), relativeVelBlockResidual, gmgSolver );
   CG->setPrintInfo( printInfo );
   auto PETScCG = std::make_shared< PETScCGSolver< P2ElementwiseAffineEpsilonOperator > >( storage, maxLevel );

   // construct pressure preconditioning operator: inverse, lumped, viscosity weighted, P1 mass matrix
   auto pPrecOp = std::make_shared< P1LumpedInvMassOperator >(
       storage,
       maxLevel,
       maxLevel
    );


   // preconditioner setup
   auto prec = std::make_shared< StokesBlockDiagonalPreconditioner< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                                                                    P1LumpedInvMassOperator> >( storage, minLevel, maxLevel, 1,  pPrecOp, LU );

   // final solver setup
   auto solver = std::make_shared< MinResSolver< P2P1ElementwiseAffineEpsilonStokesOperator > >(
       storage, minLevel, maxLevel, maxIterations, absoluteTargetResidual, prec );
   // auto solver = hyteg::MinResSolver< hyteg::P1StokesFunction< real_t >, hyteg::P1P1StokesOperator, PressurePreconditioner_T >( storage, minLevel, maxLevel, pressurePrec );
   // auto solver = hyteg::MinResSolver< hyteg::P1StokesFunction< real_t >, hyteg::P1P1StokesOperator >( storage, minLevel, maxLevel );
   solver->setPrintInfo( printInfo );
   return solver;

   //TODO smoother for MG?
   //hyteg::P1LumpedInvMassOperator massOperator( storage, minLevel, maxLevel );
   //Preconditioner_T prec( storage, minLevel, maxLevel, 2, gmgSolver );
   //DiagonalNonConstantOperator<P1ElementwiseKMassOperator, P1RowSumForm, true > ViscWeightedInvLumpedMass(storage, maxLevel, maxLevel, P1RowSumForm(forms::p1_k_mass_affine_q4));
   // BlockOperator<P2P1TaylorHoodFunction< real_t > ,P2P1TaylorHoodFunction< real_t > > stokesBlockDiagPrec(storage,minLevel,maxLevel,2,2);

   //auto CG = std::make_shared<CGSolver< P2ElementwiseAffineEpsilonOperator > >( storage, maxLevel, maxLevel);
   //CG->setPrintInfo(true);
}
*/
/// \brief Returns a geometric multigrid solver for the constant-coefficient Stokes system.
///
/// This solver performs a v-cycle and employs as smoother the inexact Uzawa method.
/// The relaxation on the A-block is performed with a forward Gauss-Seidel.
/// As coarse grid solver a pressure-preconditioned MINRES is employed.
///
/// As for all solver templates, this should not be used if performance
/// is critical. This serves as an alternative to speed up the implementation,
/// not the solution of a Stokes system solver. For fine grained control
/// the implementation can serve as a starting point.
///
/// \param storage the PrimitiveStorage that defines the domain
/// \param minLevel the coarse grid solver level
/// \param maxLevel the finest level
/// \param preSmoothingSteps number of relaxation iterations performed with the inexact Uzawa iteration
///                          before coarsening
/// \param postSmoothingSteps number of relaxation iterations performed with the inexact Uzawa iteration
///                           before after refinement
/// \param uzawaSmootherOmega relaxation parameter for the pressure preconditioner in the inexact Uzawa smoother,
///                           can be estimated by solving a related eigenvalue problem,
///                           usually 0 < omega < 1
///
template < typename StokesOperatorType >
std::shared_ptr< Solver< StokesOperatorType > > stokesGMGUzawaSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                      const uint_t&                              minLevel,
                                                                      const uint_t&                              maxLevel,
                                                                      const uint_t& preSmoothingSteps,
                                                                      const uint_t& postSmoothingSteps,
                                                                      const real_t& uzawaSmootherOmega )
{
   WALBERLA_ABORT( "Uzawa GMG solver template is not implemented for this operator." );
}

template <>
std::shared_ptr< Solver< P1P1StokesOperator > >
    stokesGMGUzawaSolver< P1P1StokesOperator >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                const uint_t&                              minLevel,
                                                const uint_t&                              maxLevel,
                                                const uint_t&                              preSmoothingSteps,
                                                const uint_t&                              postSmoothingSteps,
                                                const real_t&                              uzawaSmootherOmega )
{
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< P1P1StokesOperator, P1LumpedInvMassOperator > >(
           storage, minLevel, minLevel );
   auto pressurePreconditionedMinResSolver =
       std::make_shared< MinResSolver< P1P1StokesOperator > >( storage, minLevel, minLevel, 1000, 1e-12, pressurePreconditioner );

   auto stokesRestriction  = std::make_shared< P1P1StokesToP1P1StokesRestriction >();
   auto stokesProlongation = std::make_shared< P1P1StokesToP1P1StokesProlongation >();
   auto gaussSeidel        = std::make_shared< GaussSeidelSmoother< P1P1StokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< P1P1StokesOperator > >( storage, gaussSeidel );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< P1P1StokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaSmootherOmega );

   auto gmgSolver = std::make_shared< GeometricMultigridSolver< P1P1StokesOperator > >( storage,
                                                                                        uzawaSmoother,
                                                                                        pressurePreconditionedMinResSolver,
                                                                                        stokesRestriction,
                                                                                        stokesProlongation,
                                                                                        minLevel,
                                                                                        maxLevel,
                                                                                        preSmoothingSteps,
                                                                                        postSmoothingSteps,
                                                                                        2 );
   return gmgSolver;
}

template <>
std::shared_ptr< Solver< P2P1TaylorHoodStokesOperator > >
    stokesGMGUzawaSolver< P2P1TaylorHoodStokesOperator >( const std::shared_ptr< PrimitiveStorage >& storage,
                                                          const uint_t&                              minLevel,
                                                          const uint_t&                              maxLevel,
                                                          const uint_t&                              preSmoothingSteps,
                                                          const uint_t&                              postSmoothingSteps,
                                                          const real_t&                              uzawaSmootherOmega )
{
   auto pressurePreconditioner =
       std::make_shared< StokesPressureBlockPreconditioner< P2P1TaylorHoodStokesOperator, P1LumpedInvMassOperator > >(
           storage, minLevel, minLevel );
   auto pressurePreconditionedMinResSolver = std::make_shared< MinResSolver< P2P1TaylorHoodStokesOperator > >(
       storage, minLevel, minLevel, 1000, 1e-12, pressurePreconditioner );

   auto stokesRestriction  = std::make_shared< P2P1StokesToP2P1StokesRestriction >();
   auto stokesProlongation = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto gaussSeidel        = std::make_shared< GaussSeidelSmoother< P2P1TaylorHoodStokesOperator::VelocityOperator_T > >();
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< P2P1TaylorHoodStokesOperator > >( storage, gaussSeidel );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< P2P1TaylorHoodStokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaSmootherOmega );

   auto gmgSolver =
       std::make_shared< GeometricMultigridSolver< P2P1TaylorHoodStokesOperator > >( storage,
                                                                                     uzawaSmoother,
                                                                                     pressurePreconditionedMinResSolver,
                                                                                     stokesRestriction,
                                                                                     stokesProlongation,
                                                                                     minLevel,
                                                                                     maxLevel,
                                                                                     preSmoothingSteps,
                                                                                     postSmoothingSteps,
                                                                                     2 );
   return gmgSolver;
}

} // namespace solvertemplates
} // namespace hyteg