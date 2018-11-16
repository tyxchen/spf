# Reference: from https://www.youtube.com/watch?v=AmM0VEgfYWo

set( CMAKE_C_FLAGS "-Wall -g" )


file( GLOB SRCS *.cpp *.hpp)

add_executable( spf ${SRCS} )

find_package (PkgConfig REQUIRED)
find_package(BOOST 1.66 COMPONENTS regex  REQUIRED)
find_package (GSL REQUIRED)
target_link_libraries(spf GSL::gsl GSL::gslcblas)
