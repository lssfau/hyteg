#include <complex>
#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVersion.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/composites/P1P0StokesFunction.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/egfunctionspace/EGConvTestUtils.hpp"

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


    template<typename StokesOperatorType>
    void MultiSinker(const std::string &name,
                     const uint_t &level,
                     const uint_t &nxy,
                     const uint_t &nSinkers,
                     const real_t &DR,
                     const real_t &delta,
                     const real_t &omega) {
        using StokesFunctionType = typename StokesOperatorType::srcType;
        using StokesFunctionNumeratorType = typename StokesFunctionType::template FunctionType<idx_t>;
        real_t visc_min = std::pow(DR, -0.5);
        real_t visc_max = std::pow(DR, 0.5);

        // storage and domain
        //auto                  meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ),MeshInfo::CRISSCROSS,nxy, nxy);
        auto meshInfo = MeshInfo::meshCuboid(Point3D({0, 0, 0}), Point3D({1, 1, 1}), nxy, nxy, nxy);
        SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        hyteg::loadbalancing::roundRobin(setupStorage);
        std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, 1);

        // some general info about the testcase
        StokesFunctionNumeratorType numerator("numerator", storage, level, level);
        numerator.enumerate(level);
        uint_t globalDoFs = numberOfGlobalDoFs(numerator, level);
        WALBERLA_LOG_INFO_ON_ROOT("### Computing " << name << " ###");
        WALBERLA_LOG_INFO_ON_ROOT("Global DoFs: " << globalDoFs);
        WALBERLA_LOG_INFO_ON_ROOT("#Sinkers: " << nSinkers << ", visc_max: " << visc_max);

        // function setup
        StokesFunctionType x("x", storage, level, level);
        StokesFunctionType btmp("btmp", storage, level, level);
        StokesFunctionType b("b", storage, level, level);
        StokesFunctionType residuum("res", storage, level, level);
        std::function<real_t(const hyteg::Point3D &)> zero = [](const hyteg::Point3D &) { return real_c(0); };
        x.uvw().interpolate(zero, level);
        x.p().interpolate(zero, level);

        // generate sinker centers randomly (without them exiting the domain)
        std::uniform_real_distribution<real_t> unif(0.0, 1.0);
        std::default_random_engine re;
        re.seed(25151);
        //re.seed( 213512512 );

        std::vector<Point3D> centers;
        for (uint_t c = 0; c < nSinkers; c++) {
            centers.push_back(Point3D({unif(re), unif(re), unif(re)}));
            //centers.push_back(Point3D({unif(re), unif(re), 0}));
        }

        // characteristic function for sinkers
        std::function<real_t(const hyteg::Point3D &)> Xi = [centers, delta, omega](const hyteg::Point3D &xx) {
            real_t val = 1;
            for (auto &c: centers) {
                auto distance = c - xx;
                val *= 1 - exp(-delta * std::pow(std::max(0.0, distance.norm() - omega / 2), 2));
            }
            return val;
        };

        // viscosity function
        std::function<real_t(const hyteg::Point3D &)> viscosity = [visc_max, visc_min, Xi](const hyteg::Point3D &xx) {
            return (visc_max - visc_min) * (1 - Xi(xx)) + visc_min;
        };

        // right hand side: "force sinkers downward": negative v velocity at sinker locations
        std::function<real_t(const hyteg::Point3D &)> rhsV = [Xi](const hyteg::Point3D &xx) {
            return 10 * (Xi(xx) - 1);
        };
        btmp.uvw().interpolate({zero, zero, rhsV}, level, All);
        //btmp.uvw().interpolate({zero, rhsV}, level, All);

        if constexpr (hyteg::dg::eg::isP2P1Discr<StokesOperatorType>()) {
            P2ConstantMassOperator M_vel(storage, level, level);
            M_vel.apply(btmp.uvw()[2], b.uvw()[2], level, All);
            //M_vel.apply(btmp.uvw()[1], b.uvw()[1], level, All);

        } else if constexpr (hyteg::dg::eg::isEGP0Discr<StokesOperatorType>()) {
            hyteg::dg::eg::EGMassOperator M_vel(storage, level, level);
            M_vel.apply(btmp.uvw(), b.uvw(), level, All, Replace);
        }


        // operator setup
        StokesOperatorType Op(storage, level, level, viscosity);

        // Visualization
        VTKOutput vtkOutput("../../output", name, storage);
        vtkOutput.add(x.uvw());
        vtkOutput.add(x.p());
        vtkOutput.add(b.uvw());
        vtkOutput.add(b.p());
        if constexpr (hyteg::dg::eg::isEGP0Discr<StokesOperatorType>()) {
            vtkOutput.add(*x.uvw().getConformingPart());
            vtkOutput.add(*x.uvw().getDiscontinuousPart());
        }

        vtkOutput.write(level, 0);
        // Initial errors and residual
        //x.interpolate(zero, level, hyteg::DirichletBoundary);
        //x.interpolate([&unif, &re](const hyteg::Point3D &) { return unif(re); }, level, hyteg::Inner);
        Op.apply(x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary);
        residuum.assign({1.0, -1.0}, {b, btmp}, level, hyteg::Inner | hyteg::NeumannBoundary);
        real_t residuum_l2_1 =
                std::sqrt(residuum.dotGlobal(residuum, level, hyteg::Inner | hyteg::NeumannBoundary) /
                          (real_t) globalDoFs);
        WALBERLA_LOG_INFO_ON_ROOT("initial residual = " << residuum_l2_1);


        vtkOutput.write(level, 1);
        // Solve
        PETScLUSolver<StokesOperatorType> LU(storage, level, numerator);
        LU.solve(Op, x, b, level);

        // subtract mean
        if constexpr (hyteg::dg::eg::isEGP0Discr<StokesOperatorType>()) {
            hyteg::dg::projectMean(x.p(), level);
        } else if constexpr (hyteg::dg::eg::isP2P1Discr<StokesOperatorType>()) {
            hyteg::vertexdof::projectMean(x.p(), level);
        }

        // final residual
        Op.apply(x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary);
        residuum.assign({1.0, -1.0}, {b, btmp}, level);
        residuum_l2_1 =
                std::sqrt(residuum.dotGlobal(residuum, level, hyteg::Inner | hyteg::NeumannBoundary) /
                          (real_t) globalDoFs);
        WALBERLA_LOG_INFO_ON_ROOT("final residual = " << residuum_l2_1);

        vtkOutput.write(level, 2);
    }

} // namespace hyteg

using namespace hyteg;

int main(int argc, char *argv[]) {
    walberla::MPIManager::instance()->useWorldComm();
    PETScManager petscManager(&argc, &argv);

    MultiSinker<P2P1ElementwiseAffineEpsilonStokesOperator>("MultiSinker_P2P1", 3, 1, 4, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0", 3, 1, 4, 1000, 200, 0.1);
/*
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s1_lvl2", 2, 1, 1, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s1_lvl3", 3, 1, 1, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s1_lvl4", 4, 1, 1, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s2_lvl2", 2, 1, 2, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s2_lvl3", 3, 1, 2, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s2_lvl4", 4, 1, 2, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s4_lvl2", 2, 1, 4, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s4_lvl3", 3, 1, 4, 1000, 200, 0.1);
    MultiSinker<hyteg::dg::eg::EGP0EpsilonStokesOperator>("MultiSinker_EGP0_s4_lvl4", 4, 1, 4, 1000, 200, 0.1);
*/
    return EXIT_SUCCESS;
}
