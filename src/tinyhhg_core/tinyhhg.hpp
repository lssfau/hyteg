#ifndef TINYHHG_HPP
#define TINYHHG_HPP

#include "vtkwriter.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1CoefficientOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ElementwiseOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1PolynomialOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1BlendingOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1BlendingOperatorNew.hpp"
#include "tinyhhg_core/p1functionspace/P1PolynomialBlendingOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1HelperFunctions.hpp"
#include "p1functionspace/P1DataHandling.hpp"

#include "bubblefunctionspace/BubbleFunction.hpp"
#include "bubblefunctionspace/BubbleOperator.hpp"

#include "tinyhhg_core/dgfunctionspace/DGMemory.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFunction.hpp"
#include "tinyhhg_core/dgfunctionspace/DGUpwindOperator.hpp"

#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"

#include "types/pointnd.hpp"
#include "types/matrix.hpp"
#include "types/flags.hpp"
#include "tinyhhg_core/polynomial/Polynomial2D.hpp"
#include "tinyhhg_core/polynomial/LSQPInterpolator.hpp"

#include "solvers/cgsolver.hpp"
#include "solvers/minressolver.hpp"
#include "solvers/gmultigrid.hpp"
#include "solvers/uzawasolver.hpp"

#include "solvers/preconditioners/JacobiPreconditioner.hpp"
#include "solvers/preconditioners/GaussSeidelPreconditioner.hpp"
#include "solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"

#include "composites/p1stokesfunction.hpp"
#include "composites/p1blocklaplaceoperator.hpp"
#include "composites/p1stokesoperator.hpp"
#include "composites/P1CoefficientStokesOperator.hpp"
#include "composites/P1BlendingStokesOperator.hpp"
#include "composites/P2P1TaylorHoodFunction.hpp"
#include "composites/P2P1TaylorHoodStokesOperator.hpp"

#include "composites/P1BubbleFunctionSpace/P1BubbleFunction.hpp"
#include "composites/P1BubbleFunctionSpace/P1BubbleOperator.hpp"

#include "mixedoperators/BubbleToP1/BubbleToP1Operator.hpp"
#include "mixedoperators/P1ToBubble/P1ToBubbleOperator.hpp"

#include "composites/ministokesfunction.hpp"
#include "composites/ministokesoperator.hpp"

#include "primitivestorage/PrimitiveStorage.hpp"
#include "primitivestorage/SetupPrimitiveStorage.hpp"
#include "primitivestorage/Visualization.hpp"
#include "primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "primitivestorage/loadbalancing/DistributedBalancer.hpp"

#include "mesh/MeshInfo.hpp"

#include "communication/PackInfo.hpp"
#include "communication/DoFSpacePackInfo.hpp"
#include "communication/BufferedCommunication.hpp"

#include "petsc/PETScManager.hpp"
#include "petsc/PETScSparseMatrix.hpp"
#include "petsc/PETScVector.hpp"
#include "petsc/PETScLUSolver.hpp"
#include "petsc/PETScPreconditioner.hpp"

#include "FunctionTraits.hpp"
#include "MeshQuality.hpp"

#endif /* TINYHHG_HPP */
