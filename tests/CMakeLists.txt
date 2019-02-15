file( GLOB TESTS *.cpp *.hpp)

include_directories(${PROJECT_SOURCE_DIR}/src ${OPENMP_CXX_INCLUDE_DIRS})

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_executable(tests ${TESTS})
target_link_libraries(tests SPF OpenMP::OpenMP_CXX) # link to SPF library

add_test(SPFTests tests) # will return non-zero code if any of the tests fail
