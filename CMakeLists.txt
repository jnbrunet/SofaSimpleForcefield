cmake_minimum_required(VERSION 3.15)
project(SofaSimpleForcefield VERSION 1.0)

find_package(Sofa.Core REQUIRED)
find_package(Eigen3 REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_library(${PROJECT_NAME} SHARED init.cpp SVKElasticForcefield.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core Eigen3::Eigen)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
)