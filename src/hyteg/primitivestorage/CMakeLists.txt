target_sources( hyteg
    PRIVATE
    PrimitiveStorage.cpp     
    PrimitiveStorage.hpp
    SetupPrimitiveStorage.cpp
    SetupPrimitiveStorage.hpp
    SubStorage.hpp
    Visualization.hpp
)
add_subdirectory( loadbalancing )
