include_guard(GLOBAL)

# Add a build target which launches an OpenFPM executable with the MPI and
# runtime-library environment recorded by the installed OpenFPM package.
function(openfpm_add_mpi_run_target executable_target)
    if (NOT TARGET "${executable_target}")
        message(FATAL_ERROR
            "openfpm_add_mpi_run_target(): '${executable_target}' is not a CMake target")
    endif()

    set(_openfpm_run_one_value
        PROCESSES TARGET_NAME WORKING_DIRECTORY)
    set(_openfpm_run_multi_value
        ARGS MPI_ARGS)
    cmake_parse_arguments(OPENFPM_RUN ""
        "${_openfpm_run_one_value}"
        "${_openfpm_run_multi_value}"
        ${ARGN})

    if (OPENFPM_RUN_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "openfpm_add_mpi_run_target(): unknown arguments: ${OPENFPM_RUN_UNPARSED_ARGUMENTS}")
    endif()

    if (NOT DEFINED OPENFPM_MPI_PROCESSES)
        set(OPENFPM_MPI_PROCESSES 1 CACHE STRING
            "Number of MPI processes used by OpenFPM run targets")
    endif()
    if (NOT OPENFPM_RUN_PROCESSES)
        set(OPENFPM_RUN_PROCESSES "${OPENFPM_MPI_PROCESSES}")
    endif()
    if (NOT OPENFPM_RUN_TARGET_NAME)
        set(OPENFPM_RUN_TARGET_NAME "run_${executable_target}")
    endif()
    if (NOT OPENFPM_RUN_WORKING_DIRECTORY)
        set(OPENFPM_RUN_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    endif()

    set(_openfpm_mpiexec "${OPENFPM_MPIEXEC_EXECUTABLE}")
    if (NOT _openfpm_mpiexec OR NOT EXISTS "${_openfpm_mpiexec}")
        set(_openfpm_mpiexec "${MPIEXEC_EXECUTABLE}")
    endif()
    if (NOT _openfpm_mpiexec)
        message(FATAL_ERROR
            "openfpm_add_mpi_run_target(): no MPI launcher is available")
    endif()

    set(_openfpm_runtime_dirs)
    foreach(_openfpm_runtime_dir IN LISTS OPENFPM_RUNTIME_LIBRARY_DIRS)
        if (IS_DIRECTORY "${_openfpm_runtime_dir}")
            list(APPEND _openfpm_runtime_dirs "${_openfpm_runtime_dir}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _openfpm_runtime_dirs)

    if (WIN32)
        list(JOIN _openfpm_runtime_dirs ";" _openfpm_runtime_path)
        string(REPLACE ";" "\\;" _openfpm_runtime_path "${_openfpm_runtime_path}")
        set(_openfpm_runtime_environment
            "PATH=${_openfpm_runtime_path}\;$ENV{PATH}")
    elseif(APPLE)
        list(JOIN _openfpm_runtime_dirs ":" _openfpm_runtime_path)
        set(_openfpm_runtime_environment
            "DYLD_LIBRARY_PATH=${_openfpm_runtime_path}:$ENV{DYLD_LIBRARY_PATH}")
    else()
        list(JOIN _openfpm_runtime_dirs ":" _openfpm_runtime_path)
        set(_openfpm_runtime_environment
            "LD_LIBRARY_PATH=${_openfpm_runtime_path}:$ENV{LD_LIBRARY_PATH}")
    endif()

    add_custom_target("${OPENFPM_RUN_TARGET_NAME}"
        COMMAND "${CMAKE_COMMAND}" -E env
            "${_openfpm_runtime_environment}"
            "PURE_PYTHON=1"
            "${_openfpm_mpiexec}"
            "${MPIEXEC_NUMPROC_FLAG}" "${OPENFPM_RUN_PROCESSES}"
            ${MPIEXEC_PREFLAGS}
            ${OPENFPM_RUN_MPI_ARGS}
            "$<TARGET_FILE:${executable_target}>"
            ${MPIEXEC_POSTFLAGS}
            ${OPENFPM_RUN_ARGS}
        DEPENDS "${executable_target}"
        WORKING_DIRECTORY "${OPENFPM_RUN_WORKING_DIRECTORY}"
        USES_TERMINAL
        COMMAND_EXPAND_LISTS
        VERBATIM)
endfunction()
