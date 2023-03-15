cmake_minimum_required(VERSION 3.0.0)
project(simple_method VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories(../source_dir/build_dir)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(simple_method main.cpp vars.hpp cell.hpp f.hpp)

target_link_libraries(simple_method Eigen3::Eigen)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
