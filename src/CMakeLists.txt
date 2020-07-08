file( GLOB SRCS *.cpp *.hpp)
file(GLOB PUBLIC_INCLUDE *.hpp)
add_library(SPF STATIC ${SRCS})
set_target_properties(SPF PROPERTIES PUBLIC_HEADER ${PUBLIC_INCLUDE})
target_include_directories(SPF
    PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/include")

find_package (PkgConfig REQUIRED)
find_package(Boost 1.67 COMPONENTS regex REQUIRED)
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

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include/spf")
foreach(header ${PUBLIC_INCLUDE})
    configure_file(${header} "${CMAKE_CURRENT_BINARY_DIR}/include/spf/" COPYONLY)
endforeach()

install(TARGETS SPF
    EXPORT SPFConfig
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/spf)
export(TARGETS SPF FILE "${CMAKE_CURRENT_BINARY_DIR}/SPFConfig.cmake")
export(PACKAGE SPF)

