#ifndef TINYHHG_HPP
#define TINYHHG_HPP

#include "vtkwriter.hpp"

#include "tinyhhg_core/FunctionMemory.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "p1functionspace/P1DataHandling.hpp"

#include "bubblefunctionspace/BubbleFunction.hpp"
#include "bubblefunctionspace/BubbleOperator.hpp"

#include "types/pointnd.hpp"
#include "types/flags.hpp"

#include "solvers/cgsolver.hpp"
#include "solvers/minressolver.hpp"

#include "solvers/preconditioners/JacobiPreconditioner.hpp"
#include "solvers/preconditioners/GaussSeidelPreconditioner.hpp"

#include "composites/p1stokesfunction.hpp"
#include "composites/p1blocklaplaceoperator.hpp"
#include "composites/p1stokesoperator.hpp"

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
#include "communication/BufferedCommunication.hpp"

#include "petsc/PETScManager.hpp"
#include "petsc/PETScSparseMatrix.hpp"
#include "petsc/PETScVector.hpp"
#include "petsc/PETScLUSolver.hpp"
#include "petsc/PETScPreconditioner.hpp"



#endif /* TINYHHG_HPP */
