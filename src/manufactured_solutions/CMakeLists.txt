add_library( manufactured_solutions INTERFACE )
target_link_libraries( manufactured_solutions INTERFACE hyteg )

target_sources( manufactured_solutions
    PRIVATE
    StokesAnalyticalExpressions.hpp
)
