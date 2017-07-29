#ifndef TINYHHG_HPP
#define TINYHHG_HPP

#include "support.hpp"

#include "vtkwriter.hpp"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Memory.hpp"
#include "p1functionspace/P1DataHandling.hpp"

#include "bubblefunctionspace/BubbleFunction.hpp"
//#include "bubblefunctionspace/BubbleOperator.hpp"

//#include "mixedoperators/p1_to_p1bubble_operator.hpp"
//#include "mixedoperators/p1bubble_to_p1_operator.hpp"

#include "types/pointnd.hpp"
#include "types/flags.hpp"

#include "solvers/cgsolver.hpp"
#include "solvers/minressolver.hpp"

#include "composites/p1stokesfunction.hpp"
#include "composites/p1blocklaplaceoperator.hpp"
#include "composites/p1stokesoperator.hpp"

#include "composites/P1BubbleFunctionSpace/P1BubbleFunction.hpp"

//#include "composites/ministokesfunction.hpp"
//#include "composites/ministokesoperator.hpp"

#include "primitivestorage/PrimitiveStorage.hpp"
#include "primitivestorage/SetupPrimitiveStorage.hpp"
#include "primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "mesh/MeshInfo.hpp"

#include "communication/PackInfo.hpp"
#include "communication/BufferedCommunication.hpp"


#endif /* TINYHHG_HPP */
