# Works, but sets options for all sources of target hyteg :(
# target_compile_options( hyteg PRIVATE -Wno-shadow -Wno-unused-variable )

target_sources( hyteg
    PRIVATE
    P2PlusBubbleElementwiseDiffusion_float64.cpp
    P2PlusBubbleElementwiseDiffusion_float64.hpp
    P2PlusBubbleElementwiseMass_float64.cpp
    P2PlusBubbleElementwiseMass_float64.hpp
)
