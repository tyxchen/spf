#file( GLOB SRCS *.cpp *.hpp ${PROJECT_SOURCE_DIR}/include/*.hpp)
#include_directories(${PROJECT_SOURCE_DIR}/include)

file( GLOB SRCS *.cpp *.hpp)
add_library(SPF STATIC ${SRCS})

find_package (PkgConfig REQUIRED)
find_package(BOOST 1.66 COMPONENTS regex  REQUIRED)
find_package (GSL REQUIRED)
target_link_libraries(SPF GSL::gsl GSL::gslcblas)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

#file( GLOB PUBLIC_INCLUDE ${PROJECT_SOURCE_DIR}/include/*.hpp)
file( GLOB PUBLIC_INCLUDE *.hpp)
install(FILES ${PUBLIC_INCLUDE} DESTINATION include/spf)
install(TARGETS SPF DESTINATION bin) 
