/*
* Copyright (c) 2017-2022 Nils Kohl.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
typedef std::function<real_t(const hyteg::PointND<real_t, 3> &p)> ScalarLambda;
namespace hyteg {

    void EGApplyTest(ScalarLambda srcLambda,
                     const std::string &testName,
                     uint_t level,
                     const MeshInfo &meshInfo,
                     real_t eps,
                     bool writeVTK = false) {
        using namespace dg::eg;

        SetupPrimitiveStorage setupStorage(meshInfo,
                                           walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        auto storage = std::make_shared<PrimitiveStorage>(setupStorage, 1);

        EGFunction<real_t> src("src", storage, level, level);
        EGFunction<real_t> tmp("tmp", storage, level, level);
        EGFunction<real_t> hytegDst("hytegDst", storage, level, level);
        EGFunction<real_t> petscDst("petscDst", storage, level, level);
        EGFunction<real_t> err("error", storage, level, level);
        EGFunction<idx_t> numerator("numerator", storage, level, level);

        EGSIPGLaplaceOperator L(storage, level, level);

        // PETSc apply
        PETScVector<real_t, EGFunction> srcPetscVec;
        PETScVector<real_t, EGFunction> dstPetscVec;
        PETScSparseMatrix<EGSIPGLaplaceOperator> L_Matrix;
        PETScSparseMatrix<EGMassOperator> M_Matrix;

        std::function<real_t(const hyteg::Point3D &)> srcFunction = srcLambda;
        src.interpolate({srcFunction, srcFunction, srcFunction}, level, All);
        src.getDiscontinuousPart()->interpolate(srcFunction, level, All);

        numerator.copyBoundaryConditionFromFunction(src);
        numerator.enumerate(level);

        srcPetscVec.createVectorFromFunction(src, numerator, level);
        dstPetscVec.createVectorFromFunction(petscDst, numerator, level);
        L_Matrix.createMatrixFromOperator(L, level, numerator);
        // L_Matrix.print( "EGApplyTest_L.m", false, PETSC_VIEWER_ASCII_MATLAB );
        //  M_Matrix.print( "EGApplyTest_M.m", false, PETSC_VIEWER_ASCII_MATLAB );
        L.apply(src, hytegDst, level, All, Replace);

        // WALBERLA_CHECK( L_Matrix.isSymmetric() );

        MatMult(L_Matrix.get(), srcPetscVec.get(), dstPetscVec.get());

        dstPetscVec.createFunctionFromVector(petscDst, numerator, level);

        // compare
        err.assign({1.0, -1.0}, {hytegDst, petscDst}, level, All);
        WALBERLA_LOG_INFO_ON_ROOT("||e_disc|| = "
                                          << sqrt(err.getDiscontinuousPart()->dotGlobal(*err.getDiscontinuousPart(),
                                                                                        level, All) /
                                                  real_c(numberOfGlobalDoFs(*err.getDiscontinuousPart(), level)))
                                          << ", ||e_conf|| = "
                                          << sqrt(err.getConformingPart()->dotGlobal(*err.getConformingPart(), level,
                                                                                     All) /
                                                  real_c(numberOfGlobalDoFs(*err.getConformingPart(), level))));

        if (writeVTK) {
            VTKOutput vtk("../../output", testName, storage);

            vtk.add(src);
            vtk.add(*src.getConformingPart());
            vtk.add(*src.getDiscontinuousPart());
            vtk.add(hytegDst);
            vtk.add(*hytegDst.getConformingPart());
            vtk.add(*hytegDst.getDiscontinuousPart());
            vtk.add(petscDst);
            vtk.add(*petscDst.getConformingPart());
            vtk.add(*petscDst.getDiscontinuousPart());
            vtk.add(err);
            vtk.add(*err.getConformingPart());
            vtk.add(*err.getDiscontinuousPart());
            vtk.write(level, 0);
        }

        auto maxMag = err.getMaxMagnitude(level);

        WALBERLA_LOG_INFO_ON_ROOT(testName << ": ||e||_max = " << maxMag);

        WALBERLA_CHECK_LESS(maxMag, eps);
    }

    void EGApplyTestDiv(const std::string &testName, uint_t level, const MeshInfo &meshInfo, uint_t source,
                        bool writeVTK = false) {
        using namespace dg::eg;

        // storage declaration
        SetupPrimitiveStorage setupStorage(meshInfo,
                                           walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        auto storage = std::make_shared<PrimitiveStorage>(setupStorage, 1);

        // EG function declarations
        EGFunction<real_t> srcEG("srcEG", storage, level, level, BoundaryCondition::create0123BC());
        P0Function<real_t> hytegDst("hytegDst", storage, level, level, BoundaryCondition::createAllInnerBC());
        P0Function<real_t> petscDst("petscDst", storage, level, level, BoundaryCondition::createAllInnerBC());
        P0Function<real_t> err("error", storage, level, level, BoundaryCondition::createAllInnerBC());
        EGFunction<idx_t> numeratorP1Vec("numeratorP1Vec", storage, level, level, BoundaryCondition::create0123BC());
        P0Function<idx_t> numeratorP0("numeratorP0", storage, level, level, BoundaryCondition::createAllInnerBC());

        // P2P1, P1P0 div for comparison
        P2VectorFunction<real_t> srcP2("srcP2", storage, level, level);
        P1Function<real_t> dstP1("dstP1", storage, level, level);
        P1VectorFunction<real_t> srcP1("srcP1", storage, level, level);
        P0Function<real_t> dstP0("dstP0", storage, level, level);

        // interpolate source
        if (source == 0) {
            srcEG.interpolate(1, level, All);
            srcP2.interpolate(1, level, All);
            srcP1.interpolate(1, level, All);
        } else if (source == 1) {
            srcEG.interpolate(1, level, All);
            srcEG.getDiscontinuousPart()->interpolate(1, level, All);
            srcP2.interpolate(1, level, All);
            srcP1.interpolate(1, level, All);
        }else{
            if (storage->hasGlobalCells()) {
                srcEG.interpolate({[](const Point3D &xx) -> real_t { return -real_c(4) * std::cos(real_c(4) * xx[2]); },
                                   [](const Point3D &xx) -> real_t { return real_c(8) * std::cos(real_c(8) * xx[0]); },
                                   [](const Point3D &xx) -> real_t {
                                       return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                   }},
                                  level,
                                  All);

                srcP2.interpolate({[](const Point3D &xx) -> real_t { return -real_c(4) * std::cos(real_c(4) * xx[2]); },
                                   [](const Point3D &xx) -> real_t { return real_c(8) * std::cos(real_c(8) * xx[0]); },
                                   [](const Point3D &xx) -> real_t {
                                       return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                   }},
                                  level,
                                  All);
                srcP1.interpolate({[](const Point3D &xx) -> real_t { return -real_c(4) * std::cos(real_c(4) * xx[2]); },
                                   [](const Point3D &xx) -> real_t { return real_c(8) * std::cos(real_c(8) * xx[0]); },
                                   [](const Point3D &xx) -> real_t {
                                       return -real_c(2) * std::cos(real_c(2) * xx[1]);
                                   }},
                                  level,
                                  All);
            } else {
                srcEG.interpolate({[](const hyteg::Point3D &xx) { return real_c(20) * xx[0] * std::pow(xx[1], 3.0); },
                                   [](const hyteg::Point3D &xx) {
                                       return real_c(5) * std::pow(xx[0], 4.0) - real_c(5) * std::pow(xx[1], 4.0);
                                   }},
                                  level,
                                  All);
                srcP2.interpolate({[](const hyteg::Point3D &xx) { return real_c(20) * xx[0] * std::pow(xx[1], 3.0); },
                                   [](const hyteg::Point3D &xx) {
                                       return real_c(5) * std::pow(xx[0], 4.0) - real_c(5) * std::pow(xx[1], 4.0);
                                   }},
                                  level,
                                  All);
                srcP1.interpolate({[](const hyteg::Point3D &xx) { return real_c(20) * xx[0] * std::pow(xx[1], 3.0); },
                                   [](const hyteg::Point3D &xx) {
                                       return real_c(5) * std::pow(xx[0], 4.0) - real_c(5) * std::pow(xx[1], 4.0);
                                   }},
                                  level,
                                  All);
            }
        }

        // EGP0 and P2P1 divergence operators
        EGToP0DivOperator DivEGP0(storage, level, level);
        P2ToP1ConstantDivOperator DivP2P1(storage, level, level);
        P1ToP0ConstantDivOperator DivP1P0(storage, level, level);

        // PETSc vectors for EGP0 hyteg <-> petsc apply comparison
        PETScVector<real_t, EGFunction> srcPetscVec;
        PETScVector<real_t, P0Function> dstPetscVec;
        PETScSparseMatrix<EGToP0DivOperator> DivMatrix;


        numeratorP1Vec.enumerate(level);
        numeratorP0.enumerate(level);
        srcPetscVec.createVectorFromFunction(srcEG, numeratorP1Vec, level);
        dstPetscVec.createVectorFromFunction(petscDst, numeratorP0, level);

        // EG petsc assembly and apply
        DivMatrix.createMatrixFromOperator(DivEGP0, level, numeratorP1Vec, numeratorP0);
        MatMult(DivMatrix.get(), srcPetscVec.get(), dstPetscVec.get());
        dstPetscVec.createFunctionFromVector(petscDst, numeratorP0, level);


        // numerators for EG assembly
        //numeratorP1Vec.enumerate(level);
        //numeratorP0.enumerate(level);

        // hyteg applies
        DivEGP0.apply(srcEG, hytegDst, level, All, Replace);
        DivP2P1.apply(srcP2, dstP1, level, All, Replace);
        DivP1P0.apply(srcP1, dstP0, level, All, Replace);


        // compare EG hyteg <-> petsc apply
        err.assign({1.0, -1.0}, {hytegDst, petscDst}, level, All);

        // write for paraview
        if (writeVTK) {
            VTKOutput vtk("/mnt/c/Users/Fabia/OneDrive/Desktop/hyteg_premerge/hyteg-build/output", testName, storage);
            vtk.add(srcEG);
            vtk.add(*srcEG.getConformingPart());
            vtk.add(*srcEG.getDiscontinuousPart());
            vtk.add(hytegDst);
            vtk.add(petscDst);
            vtk.add(err);

            vtk.add(srcP2);
            vtk.add(dstP1);

            vtk.add(srcP1);
            vtk.add(dstP0);
            vtk.write(level, 0);
        }

        WALBERLA_LOG_INFO_ON_ROOT(testName << std::setprecision (15) <<": ||e|| ="
                                           << sqrt(err.dotGlobal(err, level)) / real_c(numberOfGlobalDoFs(err, level)));
        WALBERLA_LOG_INFO_ON_ROOT(testName << std::setprecision (15) <<": ||hytegDst(EGP0)|| ="
                                           << sqrt(hytegDst.dotGlobal(hytegDst, level)) /
                                              real_c(numberOfGlobalDoFs(hytegDst, level)));
        WALBERLA_LOG_INFO_ON_ROOT(testName << std::setprecision (15) <<": ||petscDst(EGP0)|| ="
                                           << sqrt(petscDst.dotGlobal(petscDst, level)) /
                                              real_c(numberOfGlobalDoFs(petscDst, level)));
        WALBERLA_LOG_INFO_ON_ROOT(testName << ": ||dst(P2P1)|| ="
                                           << sqrt(dstP1.dotGlobal(dstP1, level)) /
                                              real_c(numberOfGlobalDoFs(dstP1, level)));
        WALBERLA_LOG_INFO_ON_ROOT(testName << std::setprecision (15) <<": ||dst(P1P0)|| ="
                                           << sqrt(dstP0.dotGlobal(dstP0, level)) /
                                              real_c(numberOfGlobalDoFs(dstP0, level)));

    }

    void EGApplyTestDivt(const std::string &testName, uint_t level, const MeshInfo &meshInfo, bool writeVTK = false) {
        using namespace dg::eg;

        // storage setup
        SetupPrimitiveStorage setupStorage(meshInfo,
                                           walberla::uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
        setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);
        auto storage = std::make_shared<PrimitiveStorage>(setupStorage, 1);


        P0Function<real_t> src("src", storage, level, level);
        EGFunction<real_t> hytegDst("hytegDst", storage, level, level);
        EGFunction<real_t> petscDst("petscDst", storage, level, level);
        EGFunction<real_t> err("error", storage, level, level);
        EGFunction<idx_t> numeratorP1Vec("numeratorP1Vec", storage, level, level);
        P0Function<idx_t> numeratorP0("numeratorP0", storage, level, level);

        // P2P1, P1P0 div for comparison
        P2VectorFunction<real_t> dstP2("dstP2", storage, level, level);
        P1Function<real_t> srcP1("srcP1", storage, level, level);
        P1VectorFunction<real_t> dstP1("dstP1", storage, level, level);
        P0Function<real_t> srcP0("srcP0", storage, level, level);

        P0ToEGDivTOperator Divt(storage, level, level);
        // EGP0 and P2P1 divergence transposed operators
        P1ToP2ConstantDivTOperator DivP2P1(storage, level, level);
        //P0ToP1ConstantDivOperator DivP1P0(storage, level, level);

        // PETSc apply
        PETScVector<real_t, P0Function> srcPetscVec;
        PETScVector<real_t, EGFunction> dstPetscVec;
        PETScSparseMatrix<P0ToEGDivTOperator> Divt_Matrix;

        src.interpolate(1, level, All);

        numeratorP1Vec.enumerate(level);
        numeratorP0.enumerate(level);

        srcPetscVec.createVectorFromFunction(src, numeratorP0, level);
        dstPetscVec.createVectorFromFunction(petscDst, numeratorP1Vec, level);

        Divt_Matrix.createMatrixFromOperator(Divt, level, numeratorP0, numeratorP1Vec);
        Divt.apply(src, hytegDst, level, All, Replace);

        // WALBERLA_CHECK( L_Matrix.isSymmetric() );

        MatMult(Divt_Matrix.get(), srcPetscVec.get(), dstPetscVec.get());

        dstPetscVec.createFunctionFromVector(petscDst, numeratorP1Vec, level);

        // compare
        err.assign({1.0, -1.0}, {hytegDst, petscDst}, level, All);

        if (writeVTK) {
            VTKOutput vtk("../../output", testName, storage);

            vtk.add(src);
            vtk.add(hytegDst);
            vtk.add(petscDst);
            vtk.add(err);
            vtk.write(level, 0);
        }

        WALBERLA_LOG_INFO_ON_ROOT(testName << ": ||e||=" << sqrt(err.dotGlobal(err, level)));
        WALBERLA_LOG_INFO_ON_ROOT(testName << ": ||hytegDst|| ="
                                           << sqrt(hytegDst.dotGlobal(hytegDst, level)) /
                                              real_c(numberOfGlobalDoFs(hytegDst, level)));
        WALBERLA_LOG_INFO_ON_ROOT(testName << ": ||petscDst|| ="
                                           << sqrt(petscDst.dotGlobal(petscDst, level)) /
                                              real_c(numberOfGlobalDoFs(petscDst, level)));
    }

} // namespace hyteg

int main(int argc, char *argv[]) {
    walberla::MPIManager::instance()->initializeMPI(&argc, &argv);
    walberla::MPIManager::instance()->useWorldComm();
    hyteg::PETScManager petscManager(&argc, &argv);

    using hyteg::Point3D;
    using walberla::math::pi;
    const bool writeVTK = false;


    /*
     * // Test divT operator
    EGApplyTestDivt(
        "EGApplyDivT3D", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ), true );

    EGApplyTestDivt( "EGApplyDivT2D", 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), true );


    // test div operator
    EGApplyTestDiv("EGApplyDiv_2D_Ones", 3, hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh"), 0, true);
    EGApplyTestDiv("EGApplyDiv_2D_AllOnes", 3, hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh"), 1, true);
    EGApplyTestDiv(
            "EGApplyDiv_2D_Sinusoidal", 3,
            hyteg::MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh"), 2, true);

    EGApplyTestDiv(
            "EGApplyDiv_3D_Ones", 3,
            hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_center_at_origin_24el.msh"), 0, true);
    EGApplyTestDiv(
            "EGApplyDiv_3D_AllOnes", 3,
            hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_center_at_origin_24el.msh"), 1, true);
      EGApplyTestDiv(
            "EGApplyDiv_3D_Sinusoidal", 3,
            hyteg::MeshInfo::fromGmshFile("../../data/meshes/3D/cube_center_at_origin_24el.msh"), 2, true);
*/


    ScalarLambda srcLambda1 = []( const hyteg::Point3D& x ) { return std::sin( 3 * pi * x[0] ) * std::sin( 3 * pi * x[1] ); };
   ScalarLambda srcLambda2 = []( const hyteg::Point3D& ) { return 1; };
   ScalarLambda srcLambda3 = []( const hyteg::Point3D& x ) { return x[0] * x[0] * x[0] * std::sin( 3 * pi * x[1] ); };


   hyteg::EGApplyTest(
       srcLambda1, "tri_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "tri_2el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda1, "quad_4el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   hyteg::EGApplyTest(
       srcLambda2, "tri_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tri_2el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "quad_4el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   hyteg::EGApplyTest(
       srcLambda3, "tri_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tri_2el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/tri_2el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "quad_4el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" ), 1.0e-15, writeVTK );

   hyteg::EGApplyTest(
       srcLambda1, "tet_1el_4_src1", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda2, "tet_1el_4_src3", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   hyteg::EGApplyTest(
       srcLambda3, "tet_1el_4_src2", 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ), 1.0e-15, writeVTK );
   
   hyteg::EGApplyTest( srcLambda1,
                       "pyramid_2el_3_src1",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda2,
                       "pyramid_2el_3_src2",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda3,
                       "pyramid_2el_3_src3",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" ),
                       1.0e-15,
                       writeVTK );
   
   hyteg::EGApplyTest( srcLambda1,
                       "pyramid_4el_3_src1",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda2,
                       "pyramid_4el_3_src2",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       1.0e-15,
                       writeVTK );
   hyteg::EGApplyTest( srcLambda3,
                       "pyramid_4el_3_src3",
                       3,
                       hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_4el.msh" ),
                       1.0e-15,
                       writeVTK );

    return EXIT_SUCCESS;
}
