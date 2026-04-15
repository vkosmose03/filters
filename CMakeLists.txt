cmake_minimum_required(VERSION 3.10)
project(filters LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
    lightAI
    https://github.com/vkosmose03/lightAI.git
    GIT_TAG main
)

FetchContent_MakeAvaileble(lightAI)

add_library(filters INTERFACE)

target_include_directories(filters
    INTERFACE
    <$BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/filters/include>
    <$INSTALL_INTERFACE:include>
)

target_link_libraries(filters INTERFACE lightAI)

target_compile_features(filters INTERFACE cxx_std_17)

install(GNUInstallDirs)

install(DIRECTORY filters/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS filters EXPORT filtersTargets)

install(EXPORT filtersTargets
    FILE filtersTargets.cmake
    NAMESPACE filters::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/filters
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/filtersConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/filtersConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/filters
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/filtersConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/filtersConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/filtersConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/filters
)
