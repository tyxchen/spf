file( GLOB TESTS *.cpp *.hpp)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(tests ${TESTS})
target_link_libraries(tests SPF)
