cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)
project(repetitiveness)

# include CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

# SDSL
find_package(SDSL REQUIRED)

# set C++ build flags
set(CXX_STANDARD c++20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=gnu++23 ${GCC_WARNINGS} ${OpenMP_CXX_FLAGS} -Wstringop-overflow=0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -march=native")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

if(${BENCHMARK})
    add_compile_definitions(BENCHMARK)
endif()

add_executable(repetitiveness src/main.cpp)
target_include_directories(repetitiveness PRIVATE ${SDSL_INCLUDE_DIRS})
target_link_libraries(repetitiveness ${SDSL_LIBRARIES})
