#include <core/timing/Timer.h>

#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#ifndef HHG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHHG_BUILD_WITH_PETSC=ON"
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main(int argc, char *argv[]) {
    walberla::Environment walberlaEnv(argc, argv);
    walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::PROGRESS);
    walberla::MPIManager::instance()->useWorldComm();

    PETScManager petscManager;

    std::string meshFileName = "../../data/meshes/unitsquare_with_circular_hole.msh";

    MeshInfo meshInfo = MeshInfo::fromGmshFile(meshFileName);
    SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));

    hhg::loadbalancing::roundRobin(setupStorage);

    const size_t level = 4;

    std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

    hhg::P1Function<real_t> x("x", storage, level, level);
    hhg::P1Function<real_t> Ax("Ax", storage, level, level);
    hhg::P1Function<real_t> b("b", storage, level, level);
    hhg::P1Function<real_t> x_exact("x_exact", storage, level, level);
    hhg::P1Function<real_t> err("err", storage, level, level);

    hhg::P1ConstantLaplaceOperator A(storage, level, level);

    std::function<real_t(const hhg::Point3D &)> exact = [](const hhg::Point3D &xx) {
        return xx[0] * xx[0] - xx[1] * xx[1] + 10;
    };
    std::function<real_t(const hhg::Point3D &)> rhs = [](const hhg::Point3D &) { return 0.0; };
    std::function<real_t(const hhg::Point3D &)> ones = [](const hhg::Point3D &) { return 1.0; };

    b.interpolate(exact, level, hhg::DirichletBoundary);
    x.interpolate(exact, level, hhg::DirichletBoundary);
    x_exact.interpolate(exact, level);

    uint_t globalDoFs = hhg::numberOfGlobalDoFs<P1FunctionTag>(*storage, level);
    WALBERLA_LOG_INFO_ON_ROOT("Num dofs = " << uint_c(globalDoFs))

    PETScLUSolver<hhg::P1ConstantLaplaceOperator> solver(storage, level);

    // Petsc output
    uint_t localDoFs = hhg::numberOfLocalDoFs<P1FunctionTag>(*storage, level);
    P1Function<PetscInt> numerator("numerator", storage, level, level);
    numerator.enumerate(level);
    PETScVector<real_t, P1Function> bVec(localDoFs);
    PETScVector<real_t, P1Function> xVec(localDoFs);
    PETScVector<real_t, P1Function> exaPetscVec(localDoFs);
    PETScSparseMatrix<P1ConstantLaplaceOperator, P1Function> petscMatrix(localDoFs, globalDoFs);
    bVec.createVectorFromFunction(b, numerator, level, All);
    exaPetscVec.createVectorFromFunction(x_exact, numerator, level, All);
    petscMatrix.createMatrixFromFunction(A, level, numerator, All);
//    bVec.createVectorFromFunction(b, numerator, level, Inner);
//    exaPetscVec.createVectorFromFunction(x_exact, numerator, level, Inner);
//    petscMatrix.createMatrixFromFunction(A, level, numerator, Inner);
    petscMatrix.applyDirichletBC(numerator, level);
    bVec.print("bLU.m");
    exaPetscVec.print("exLU.m");
    petscMatrix.print("ALU.m");
    WALBERLA_LOG_INFO_ON_ROOT("Solving System")
    walberla::WcTimer timer;
    solver.solve(A, x, b, level);
    timer.end();

    xVec.createVectorFromFunction(x, numerator, level, All);
//    xVec.createVectorFromFunction(x, numerator, level, Inner);
    xVec.print("xLU.m");

    // Residual computation
    A.apply(x, Ax, level, hhg::Inner);
    err.assign({1.0, -1.0}, {b, Ax}, level, hhg::Inner);
    real_t abs_res_old = std::sqrt(err.dotGlobal(err, level, hhg::Inner));
    WALBERLA_LOG_INFO_ON_ROOT("discrete L2 residual = " << abs_res_old);


    // Petsc output DONE


    WALBERLA_LOG_INFO_ON_ROOT("time was: " << timer.last());
    err.assign({1.0, -1.0}, {x, x_exact}, level);

    real_t discr_l2_err = std::sqrt(err.dotGlobal(err, level) / (real_t) globalDoFs);

    WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);
//   WALBERLA_CHECK_LESS( discr_l2_err, 1e-14 );

    //  WALBERLA_LOG_INFO_ON_ROOT("Printing Solution")
    //  hhg::VTKWriter< P1Function >({ x, x_exact, &err }, level, "../output", "exact_solver");

    return EXIT_SUCCESS;
}
