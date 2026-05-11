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
option(EPCSAFT_ENABLE_CERES "Enable package-wide Ceres support for native thermodynamic and regression residual backends" OFF)
option(EPCSAFT_USE_SYSTEM_CERES "Use an installed Ceres Solver package instead of FetchContent" OFF)
option(EPCSAFT_ENABLE_CPPAD "Enable package-wide CppAD support for native derivative/automatic-differentiation workflows" OFF)
option(EPCSAFT_USE_SYSTEM_CPPAD "Use an installed CppAD package instead of bundled/fallback sources" OFF)
option(EPCSAFT_WARNINGS_AS_ERRORS "Treat native warnings as errors" OFF)
option(EPCSAFT_ENABLE_SANITIZERS "Enable ASAN/UBSAN in debug builds" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import includeigen; print(includeigen.get_include())"
    OUTPUT_VARIABLE EPCSAFT_EIGEN_INCLUDE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)
set(EIGEN3_INCLUDE_DIR "${EPCSAFT_EIGEN_INCLUDE}" CACHE PATH "Eigen include directory for native dependencies" FORCE)

if(EPCSAFT_ENABLE_CERES)
    if(EPCSAFT_USE_SYSTEM_CERES)
        find_package(Ceres CONFIG REQUIRED)
    else()
        include(FetchContent)
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
        set(MINIGLOG ON CACHE BOOL "" FORCE)
        FetchContent_Declare(
            ceres-solver
            GIT_REPOSITORY https://ceres-solver.googlesource.com/ceres-solver
            GIT_TAG 2.2.0
        )
        FetchContent_MakeAvailable(ceres-solver)
        if(MINGW AND TARGET ceres_internal)
            target_compile_options(ceres_internal PRIVATE -include cstdint)
        endif()
    endif()
    set(EPCSAFT_CERES_AVAILABLE ON CACHE BOOL "Ceres support is available for ePC-SAFT native code" FORCE)
else()
    set(EPCSAFT_CERES_AVAILABLE OFF CACHE BOOL "Ceres support is available for ePC-SAFT native code" FORCE)
endif()

if(EPCSAFT_ENABLE_CPPAD)
    if(EPCSAFT_USE_SYSTEM_CPPAD)
        find_path(CPPAD_INCLUDE_DIR
            NAMES cppad/cppad.hpp
            REQUIRED
        )
        set(EPCSAFT_CPPAD_INCLUDE_DIR "${CPPAD_INCLUDE_DIR}")
    else()
        include(FetchContent)
        FetchContent_Declare(
            cppad
            GIT_REPOSITORY https://github.com/coin-or/CppAD.git
            GIT_TAG 20250000.3
            SOURCE_SUBDIR cmake-do-not-add
        )
        FetchContent_MakeAvailable(cppad)
        set(compiler_has_conversion_warn 0)
        set(cppad_debug_and_release_01 0)
        set(use_cplusplus_2017_ok 1)
        set(cppad_version "20250000.3")
        set(cppad_has_adolc 0)
        set(cppad_has_colpack 0)
        set(cppad_has_eigen 1)
        set(cppad_has_ipopt 0)
        set(cppad_deprecated_01 0)
        set(cppad_boostvector 0)
        set(cppad_cppadvector 1)
        set(cppad_stdvector 0)
        set(cppad_eigenvector 0)
        set(cppad_has_gettimeofday 0)
        set(cppad_tape_addr_type "unsigned int")
        set(cppad_is_same_tape_addr_type_size_t 0)
        set(cppad_tape_id_type "unsigned int")
        set(cppad_max_num_threads 48)
        set(cppad_has_mkstemp 0)
        set(cppad_has_tmpnam_s 0)
        set(cppad_c_compiler_cmd "${CMAKE_C_COMPILER}")
        set(cppad_c_compiler_gnu_flags 0)
        set(cppad_c_compiler_msvc_flags 0)
        set(cppad_is_same_unsigned_int_size_t 0)
        set(cppad_padding_block_t "char padding_[8];")
        configure_file(
            "${cppad_SOURCE_DIR}/include/cppad/configure.hpp.in"
            "${CMAKE_CURRENT_BINARY_DIR}/cppad-config/include/cppad/configure.hpp"
            @ONLY
        )
        set(EPCSAFT_CPPAD_INCLUDE_DIR
            "${CMAKE_CURRENT_BINARY_DIR}/cppad-config/include"
            "${cppad_SOURCE_DIR}/include"
        )
        set(EPCSAFT_CPPAD_SOURCES "${cppad_SOURCE_DIR}/cppad_lib/temp_file.cpp")
    endif()
    set(EPCSAFT_CPPAD_AVAILABLE ON CACHE BOOL "CppAD headers are available for ePC-SAFT native code" FORCE)
else()
    set(EPCSAFT_CPPAD_AVAILABLE OFF CACHE BOOL "CppAD headers are available for ePC-SAFT native code" FORCE)
endif()

file(GLOB EPCSAFT_NATIVE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/autodiff/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/regression/*.cpp"
)

add_library(epcsaft_native STATIC ${EPCSAFT_NATIVE_SOURCES})
if(EPCSAFT_CPPAD_SOURCES)
    target_sources(epcsaft_native PRIVATE ${EPCSAFT_CPPAD_SOURCES})
endif()
target_include_directories(epcsaft_native
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/autodiff"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/regression"
        "${EPCSAFT_EIGEN_INCLUDE}"
)
target_compile_features(epcsaft_native PUBLIC cxx_std_17)
set_target_properties(epcsaft_native PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(EPCSAFT_ENABLE_CERES)
    target_link_libraries(epcsaft_native PUBLIC Ceres::ceres)
    target_compile_definitions(epcsaft_native PUBLIC EPCSAFT_HAS_CERES=1)
endif()
if(EPCSAFT_ENABLE_CPPAD)
    target_include_directories(epcsaft_native PUBLIC ${EPCSAFT_CPPAD_INCLUDE_DIR})
    target_compile_definitions(epcsaft_native PUBLIC EPCSAFT_HAS_CPPAD=1)
endif()

pybind11_add_module(_core src/epcsaft/bindings.cpp)
target_link_libraries(_core PRIVATE epcsaft_native)
target_include_directories(_core
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/autodiff"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/contributions"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/equilibrium"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/epcsaft/native/regression"
        "${EPCSAFT_EIGEN_INCLUDE}"
)
if(EPCSAFT_ENABLE_CPPAD)
    target_include_directories(_core PRIVATE ${EPCSAFT_CPPAD_INCLUDE_DIR})
endif()
target_compile_features(_core PRIVATE cxx_std_17)

set(EPCSAFT_NATIVE_TARGETS epcsaft_native _core)

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
