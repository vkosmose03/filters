cmake_minimum_required(VERSION 3.10)
project(filters)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/filters/include)

# Create example executable - добавляем только example.cpp
# Остальные файлы не нужны, так как их реализация в заголовках
add_executable(example
    example/example.cpp
)

# Enable warnings and debugging
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(example PRIVATE -Wall -Wextra -g)
endif()