if(NOT DEFINED EIGEN3_INCLUDE_DIR OR EIGEN3_INCLUDE_DIR STREQUAL "")
    if(DEFINED Python_EXECUTABLE)
        execute_process(
            COMMAND "${Python_EXECUTABLE}" -c "import includeigen; print(includeigen.get_include())"
            OUTPUT_VARIABLE EIGEN3_INCLUDE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            COMMAND_ERROR_IS_FATAL ANY
        )
    endif()
endif()

if(DEFINED EIGEN3_INCLUDE_DIR AND EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Core")
    set(Eigen3_FOUND TRUE)
    set(EIGEN3_FOUND TRUE)
    set(Eigen3_VERSION "3.4.0")
    set(EIGEN3_VERSION_STRING "3.4.0")
    if(NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    endif()
else()
    set(Eigen3_FOUND FALSE)
    set(EIGEN3_FOUND FALSE)
endif()
