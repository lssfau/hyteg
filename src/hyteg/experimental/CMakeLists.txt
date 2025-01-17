target_sources( hyteg
    PRIVATE
    P2PlusBubbleFunction.cpp
    P2PlusBubbleFunction.hpp
    SphericalElementFunction.hpp
)

add_subdirectory( P2PlusBubbleOperators )
set_source_files_properties( DIRECTORY P2PlusBubbleOperators PROPERTIES COMPILE_OPTIONS "-Wno-unused-variable;-Wno-shadow" )
