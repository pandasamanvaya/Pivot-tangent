cmake_minimum_required(VERSION 3.12)
project(InverseSqrtComparision VERSION 1.0)

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})
find_package(SEAL)
set (CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
    matplotlib
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        f23347fca25219d1c42cbb91608b5556814bf572
)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
    FetchContent_Populate(matplotlib)
endif()
include_directories(SYSTEM ${matplotlib_SOURCE_DIR})

add_executable(comp inv_sqrt_comp.cpp)
target_link_libraries(comp PRIVATE SEAL::seal)
target_link_libraries(comp
    PRIVATE
        ${PYTHON_LIBRARIES}
        Python3::NumPy
)