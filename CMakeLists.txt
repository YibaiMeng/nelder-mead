cmake_minimum_required(VERSION 3.5)

project(nelder-mead-test LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)

file(GLOB SOURCES
    ${CMAKE_SOURCE_DIR}/*.h 
    ${CMAKE_SOURCE_DIR}/*.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCES})

target_compile_options(${PROJECT_NAME} PUBLIC 
    -Wall -Wextra
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g3 -O0 -fno-inline>
    $<$<CONFIG:Profile>:-O3 -pg>
    )
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_EXTENSIONS OFF)
set_property(TARGET ${PROJECT_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}
    )

target_link_libraries(${PROJECT_NAME})
