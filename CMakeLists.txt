cmake_minimum_required(VERSION 3.14)
project(filters VERSION 1.0.0 LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
    lightAI
    GIT_REPOSITORY https://github.com/vkosmose03/lightAI.git
    GIT_TAG main
)

set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(lightAI)

add_library(filters INTERFACE)

target_include_directories(filters INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(filters INTERFACE lightAI::core)

target_compile_features(filters INTERFACE cxx_std_17)
