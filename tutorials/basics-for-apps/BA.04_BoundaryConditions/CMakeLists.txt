waLBerla_add_executable( NAME BA.04_BoundaryConditions
        FILES BA.04_BoundaryConditions.cpp
        DEPENDS hyteg core )
waLBerla_execute_test( NO_MODULE_LABEL NAME BA.04_BoundaryConditions )