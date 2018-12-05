file( GLOB TESTS *.cpp *.hpp)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(tests ${TESTS})
target_link_libraries(tests SPF) # link to SPF library

add_test(SPFTests tests) # will return non-zero code if any of the tests fail
