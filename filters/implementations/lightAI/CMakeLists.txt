cmake_minimum_required(VERSION 3.16)
project(lightAI CXX)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
add_library(lightAI_core
    core/layer.cpp
    core/lightAI.cpp
    core/filterMLP.cpp
)
target_include_directories(lightAI_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
)
target_link_libraries(lightAI_core PUBLIC Eigen3::Eigen)
add_library(lightAI::core ALIAS lightAI_core)
