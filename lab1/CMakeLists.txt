add_library(modular_arithmetics modular.h modular.cpp)
add_library(polynomial_arithmetics polynomial.h polynomial.cpp)

add_executable(lab1 main.cpp)

target_link_libraries(lab1 PRIVATE modular_arithmetics polynomial_arithmetics)

set_target_properties(lab1 modular_arithmetics polynomial_arithmetics
PROPERTIES
RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug/lab1/bin
ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug/lab1/lib
)