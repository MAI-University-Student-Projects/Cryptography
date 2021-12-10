# add_library(modular_arithmetics modular.hpp modular.cpp) #for header only INTERFACE key-word must be used
# add_library(polynomial_arithmetics polynomial.hpp polynomial.cpp)

add_executable(lab1 main.cpp primes.hpp modular.hpp polynomial.hpp)

# target_link_libraries(lab1 PRIVATE modular_arithmetics polynomial_arithmetics)

# set_target_properties(lab1 modular_arithmetics polynomial_arithmetics 
# PROPERTIES
# RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug/lab1/bin
# ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug/lab1/lib
# )

set_target_properties(lab1 PROPERTIES
RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug/lab1/bin
ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug/lab1/lib
)
