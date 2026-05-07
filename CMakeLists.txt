cmake_minimum_required(VERSION 3.20...4.2)

if(NOT DEFINED SKBUILD_PROJECT_NAME)
    set(SKBUILD_PROJECT_NAME epcsaft)
endif()
if(NOT DEFINED SKBUILD_PROJECT_VERSION)
    set(SKBUILD_PROJECT_VERSION 1.5.0)
endif()

project(${SKBUILD_PROJECT_NAME}
    VERSION ${SKBUILD_PROJECT_VERSION}
    LANGUAGES CXX
)

option(EPCSAFT_DEV_INPLACE "Copy the compiled extension into src/epcsaft for local development" OFF)
option(EPCSAFT_ENABLE_CERES "Enable Ceres Solver support for native equilibrium nonlinear solves" OFF)
option(EPCSAFT_USE_SYSTEM_CERES "Use an installed Ceres Solver package instead of FetchContent" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)

if(EPCSAFT_ENABLE_CERES)
    if(EPCSAFT_USE_SYSTEM_CERES)
        find_package(Ceres CONFIG REQUIRED)
    else()
        include(FetchContent)
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(
            ceres-solver
            GIT_REPOSITORY https://ceres-solver.googlesource.com/ceres-solver
            GIT_TAG 2.2.0
        )
        FetchContent_MakeAvailable(ceres-solver)
    endif()
endif()

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import includeigen; print(includeigen.get_include())"
    OUTPUT_VARIABLE EPCSAFT_EIGEN_INCLUDE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)

file(GLOB EPCSAFT_NATIVE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium/*.cpp"
)

add_library(epcsaft_native STATIC ${EPCSAFT_NATIVE_SOURCES})
target_include_directories(epcsaft_native
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium"
        "${EPCSAFT_EIGEN_INCLUDE}"
)
target_compile_features(epcsaft_native PUBLIC cxx_std_17)
set_target_properties(epcsaft_native PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(EPCSAFT_ENABLE_CERES)
    target_link_libraries(epcsaft_native PUBLIC Ceres::ceres)
    target_compile_definitions(epcsaft_native PUBLIC EPCSAFT_HAS_CERES=1)
endif()

pybind11_add_module(_core src/epcsaft/bindings.cpp)
target_link_libraries(_core PRIVATE epcsaft_native)
target_include_directories(_core
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium"
        "${EPCSAFT_EIGEN_INCLUDE}"
)
target_compile_features(_core PRIVATE cxx_std_17)

if(MSVC)
    target_compile_options(epcsaft_native PRIVATE /wd4551)
    target_compile_options(_core PRIVATE /wd4551)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(epcsaft_native PRIVATE -Wa,-mbig-obj)
    target_compile_options(_core PRIVATE -Wa,-mbig-obj)
    target_link_options(_core PRIVATE -static -static-libgcc -static-libstdc++)
endif()

install(TARGETS _core
    LIBRARY DESTINATION epcsaft
    RUNTIME DESTINATION epcsaft
)

if(EPCSAFT_DEV_INPLACE)
    add_custom_command(TARGET _core POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "$<TARGET_FILE:_core>"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/$<TARGET_FILE_NAME:_core>"
        VERBATIM
    )
endif()
