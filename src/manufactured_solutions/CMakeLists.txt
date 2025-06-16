add_library( manufactured_solutions )
target_link_libraries( manufactured_solutions PUBLIC hyteg )

target_sources( manufactured_solutions
    PRIVATE
    StokesAnalyticalExpressions.hpp
)
