file( GLOB SRCS *.cpp *.hpp)
add_library(SPF STATIC ${SRCS})

find_package (PkgConfig REQUIRED)
find_package(BOOST 1.66 COMPONENTS regex  REQUIRED)
find_package (GSL REQUIRED)
target_link_libraries(SPF GSL::gsl GSL::gslcblas)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

file( GLOB PUBLIC_INCLUDE *.hpp) # grab all hpp files
install(FILES ${PUBLIC_INCLUDE} DESTINATION include/spf) # copy it to include/spf
install(TARGETS SPF DESTINATION bin) # install .a file to bin/
