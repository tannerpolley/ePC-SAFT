cmake_minimum_required(VERSION 3.20...4.2)

if(NOT DEFINED SKBUILD_PROJECT_NAME)
    set(SKBUILD_PROJECT_NAME epcsaft)
endif()
if(NOT DEFINED SKBUILD_PROJECT_VERSION)
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml" EPCSAFT_PYPROJECT_TOML)
    string(REGEX MATCH "version[ \t]*=[ \t]*\"([^\"]+)\"" EPCSAFT_VERSION_MATCH "${EPCSAFT_PYPROJECT_TOML}")
    if(NOT CMAKE_MATCH_1)
        message(FATAL_ERROR "Could not derive SKBUILD_PROJECT_VERSION from pyproject.toml")
    endif()
    set(SKBUILD_PROJECT_VERSION "${CMAKE_MATCH_1}")
endif()

project(${SKBUILD_PROJECT_NAME}
    VERSION ${SKBUILD_PROJECT_VERSION}
    LANGUAGES CXX
)

option(EPCSAFT_DEV_INPLACE "Copy the compiled extension into src/epcsaft for local development" OFF)
option(EPCSAFT_ENABLE_CERES "Enable Ceres Solver support for native equilibrium nonlinear solves" ON)
option(EPCSAFT_USE_SYSTEM_CERES "Use an installed Ceres Solver package instead of FetchContent" OFF)
option(EPCSAFT_ENABLE_CPPAD "Enable package-wide CppAD support" ON)
option(EPCSAFT_USE_SYSTEM_CPPAD "Use installed CppAD instead of FetchContent" OFF)
option(EPCSAFT_WARNINGS_AS_ERRORS "Treat native warnings as errors" OFF)
option(EPCSAFT_ENABLE_SANITIZERS "Enable ASAN/UBSAN in debug builds" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import includeigen; print(includeigen.get_include())"
    OUTPUT_VARIABLE EPCSAFT_EIGEN_INCLUDE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)
file(TO_CMAKE_PATH "${EPCSAFT_EIGEN_INCLUDE}" EPCSAFT_EIGEN_INCLUDE_CMAKE)

if(EPCSAFT_ENABLE_CERES AND NOT Eigen3_DIR)
    set(EPCSAFT_EIGEN3_CONFIG_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated_eigen3_config")
    file(MAKE_DIRECTORY "${EPCSAFT_EIGEN3_CONFIG_DIR}")
    file(WRITE "${EPCSAFT_EIGEN3_CONFIG_DIR}/Eigen3Config.cmake"
"if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES \"${EPCSAFT_EIGEN_INCLUDE_CMAKE}\")
endif()
set(Eigen3_VERSION 3.4.0)
set(Eigen3_FOUND TRUE)
")
    file(WRITE "${EPCSAFT_EIGEN3_CONFIG_DIR}/Eigen3ConfigVersion.cmake"
"set(PACKAGE_VERSION 3.4.0)
set(PACKAGE_VERSION_COMPATIBLE TRUE)
set(PACKAGE_VERSION_EXACT FALSE)
")
    set(Eigen3_DIR "${EPCSAFT_EIGEN3_CONFIG_DIR}" CACHE PATH "Generated Eigen3 package config for bundled includeigen")
endif()

if(EPCSAFT_ENABLE_CERES)
    if(EPCSAFT_USE_SYSTEM_CERES)
        find_package(Ceres CONFIG REQUIRED)
    else()
        include(FetchContent)
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
        set(MINIGLOG ON CACHE BOOL "" FORCE)
        set(GFLAGS OFF CACHE BOOL "" FORCE)
        set(SUITESPARSE OFF CACHE BOOL "" FORCE)
        set(EIGENSPARSE OFF CACHE BOOL "" FORCE)
        set(LAPACK OFF CACHE BOOL "" FORCE)
        set(USE_CUDA OFF CACHE BOOL "" FORCE)
        set(SCHUR_SPECIALIZATIONS OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(
            ceres-solver
            GIT_REPOSITORY https://ceres-solver.googlesource.com/ceres-solver
            GIT_TAG 2.2.0
        )
        FetchContent_MakeAvailable(ceres-solver)
        if(TARGET ceres_internal AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            target_compile_options(ceres_internal PRIVATE -include cstdint)
        endif()
    endif()
endif()

if(EPCSAFT_USE_SYSTEM_CPPAD)
    set(EPCSAFT_ENABLE_CPPAD ON CACHE BOOL "Enable package-wide CppAD support" FORCE)
endif()

set(EPCSAFT_CPPAD_STATUS "disabled")
set(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR "")
set(EPCSAFT_CPPAD_SOURCE_FILE "")
set(EPCSAFT_CPPAD_LIBRARY "")
if(EPCSAFT_ENABLE_CPPAD)
    if(EPCSAFT_USE_SYSTEM_CPPAD)
        find_path(EPCSAFT_CPPAD_INCLUDE_DIR NAMES cppad/cppad.hpp REQUIRED)
        find_library(EPCSAFT_CPPAD_LIBRARY NAMES cppad_lib cppad)
    else()
        include(FetchContent)
        FetchContent_Declare(
            cppad
            GIT_REPOSITORY https://github.com/coin-or/CppAD.git
            GIT_TAG 20250000.3
            GIT_SHALLOW TRUE
            SOURCE_SUBDIR cmake-do-not-configure
        )
        FetchContent_MakeAvailable(cppad)
        set(EPCSAFT_CPPAD_INCLUDE_DIR "${cppad_SOURCE_DIR}/include")
        set(EPCSAFT_CPPAD_SOURCE_FILE "${cppad_SOURCE_DIR}/cppad_lib/temp_file.cpp")
    endif()
    if(NOT EXISTS "${EPCSAFT_CPPAD_INCLUDE_DIR}/cppad/cppad.hpp")
        message(FATAL_ERROR "EPCSAFT_ENABLE_CPPAD=ON but cppad/cppad.hpp was not found")
    endif()
    if(NOT EXISTS "${EPCSAFT_CPPAD_INCLUDE_DIR}/cppad/configure.hpp")
        if(NOT EXISTS "${EPCSAFT_CPPAD_INCLUDE_DIR}/cppad/configure.hpp.in")
            message(FATAL_ERROR "EPCSAFT_ENABLE_CPPAD=ON but cppad/configure.hpp is missing and no template was found")
        endif()
        set(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated_cppad_include")
        file(MAKE_DIRECTORY "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}/cppad")
        set(compiler_has_conversion_warn 0)
        set(cppad_debug_and_release_01 1)
        set(use_cplusplus_2017_ok 1)
        set(cppad_version "20250000.3")
        set(cppad_has_adolc 0)
        set(cppad_has_colpack 0)
        set(cppad_has_eigen 0)
        set(cppad_has_ipopt 0)
        set(cppad_deprecated_01 0)
        set(cppad_boostvector 0)
        set(cppad_cppadvector 0)
        set(cppad_stdvector 1)
        set(cppad_eigenvector 0)
        set(cppad_has_gettimeofday 0)
        set(cppad_tape_addr_type "unsigned int")
        set(cppad_is_same_tape_addr_type_size_t 0)
        set(cppad_tape_id_type "unsigned int")
        set(cppad_max_num_threads 48)
        set(cppad_has_mkstemp 0)
        set(cppad_has_tmpnam_s 0)
        set(cppad_c_compiler_cmd "${CMAKE_CXX_COMPILER}")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            set(cppad_c_compiler_gnu_flags 1)
        else()
            set(cppad_c_compiler_gnu_flags 0)
        endif()
        if(MSVC)
            set(cppad_c_compiler_msvc_flags 1)
        else()
            set(cppad_c_compiler_msvc_flags 0)
        endif()
        set(cppad_is_same_unsigned_int_size_t 0)
        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            set(cppad_padding_block_t "")
        else()
            set(cppad_padding_block_t "void* note_used_;")
        endif()
        configure_file(
            "${EPCSAFT_CPPAD_INCLUDE_DIR}/cppad/configure.hpp.in"
            "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}/cppad/configure.hpp"
            @ONLY
        )
    endif()
    set(EPCSAFT_CPPAD_STATUS "enabled_available")
endif()

file(GLOB EPCSAFT_NATIVE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/autodiff/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium/*.cpp"
)
if(EPCSAFT_CPPAD_SOURCE_FILE)
    list(APPEND EPCSAFT_NATIVE_SOURCES "${EPCSAFT_CPPAD_SOURCE_FILE}")
endif()

add_library(epcsaft_native OBJECT ${EPCSAFT_NATIVE_SOURCES})
target_include_directories(epcsaft_native
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/autodiff"
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
if(EPCSAFT_CPPAD_LIBRARY)
    target_link_libraries(epcsaft_native PUBLIC "${EPCSAFT_CPPAD_LIBRARY}")
endif()

pybind11_add_module(_core src/epcsaft/bindings.cpp)
target_link_libraries(_core PRIVATE epcsaft_native)
target_include_directories(_core
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/autodiff"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium"
        "${EPCSAFT_EIGEN_INCLUDE}"
)
target_compile_features(_core PRIVATE cxx_std_17)

set(EPCSAFT_NATIVE_TARGETS epcsaft_native _core)

foreach(target_name IN LISTS EPCSAFT_NATIVE_TARGETS)
    target_compile_definitions(${target_name} PRIVATE EPCSAFT_CPPAD_STATUS="${EPCSAFT_CPPAD_STATUS}")
    if(EPCSAFT_ENABLE_CPPAD)
        if(EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR)
            target_include_directories(${target_name} SYSTEM PRIVATE "${EPCSAFT_CPPAD_GENERATED_INCLUDE_DIR}")
        endif()
        target_include_directories(${target_name} SYSTEM PRIVATE "${EPCSAFT_CPPAD_INCLUDE_DIR}")
        target_compile_definitions(${target_name} PRIVATE EPCSAFT_HAS_CPPAD=1)
    endif()
endforeach()

foreach(target_name IN LISTS EPCSAFT_NATIVE_TARGETS)
    if(MSVC)
        target_compile_options(${target_name} PRIVATE /W4 /wd4551)
        if(EPCSAFT_WARNINGS_AS_ERRORS)
            target_compile_options(${target_name} PRIVATE /WX)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic)
        if(EPCSAFT_WARNINGS_AS_ERRORS)
            target_compile_options(${target_name} PRIVATE -Werror)
        endif()
        if(EPCSAFT_ENABLE_SANITIZERS AND CMAKE_BUILD_TYPE MATCHES "Debug")
            target_compile_options(${target_name} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
            target_link_options(${target_name} PRIVATE -fsanitize=address,undefined)
        endif()
    endif()
endforeach()

if(EPCSAFT_ENABLE_SANITIZERS AND MSVC)
    message(WARNING "EPCSAFT_ENABLE_SANITIZERS is currently only configured for GNU/Clang debug builds.")
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
