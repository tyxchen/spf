# define the project name: spf-lib
project(spf-lib)

cmake_minimum_required (VERSION 3.6)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set( CMAKE_C_FLAGS "-Wall -g" )

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)

