cmake_minimum_required(VERSION 3.16)
project(NWQEC VERSION 1.0 LANGUAGES CXX)

# =============================================================================
# Project Configuration
# =============================================================================

# Enable C++17 features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Display compiler information
message(STATUS "Using C++ compiler: ${CMAKE_CXX_COMPILER}")

# Options
option(NWQEC_ENABLE_LTO "Enable IPO/LTO in Release" ON)
option(NWQEC_ENABLE_NATIVE "Enable -march=native in Release (GCC/Clang)" OFF)
option(NWQEC_BUILD_PYTHON "Build Python bindings (pybind11)" ON)

option(NWQEC_ALLOW_NO_GMP "Allow building without GMP/MPFR (use Python gridsynth fallback)" OFF)

# =============================================================================
# Compiler Options
# =============================================================================

# Common compiler options
set(COMMON_COMPILE_OPTIONS)
if(MSVC)
    list(APPEND COMMON_COMPILE_OPTIONS /W4)
else()
    list(APPEND COMMON_COMPILE_OPTIONS -Wall -Wextra -pedantic)
endif()

# =============================================================================
# Library Targets
# =============================================================================

# Public interface for the project headers
add_library(nwqec INTERFACE)
target_include_directories(nwqec INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(nwqec INTERFACE cxx_std_17)

# =============================================================================
# Gridsynth dependency (GMP/MPFR) and target
# =============================================================================

# Try find_package first; if unavailable, create imported targets from paths
find_package(GMP QUIET)
find_package(MPFR QUIET)

if(NOT TARGET GMP::gmp OR NOT TARGET GMP::gmpxx)
    find_path(GMP_INCLUDE_DIR NAMES gmp.h)
    find_library(GMP_LIBRARY NAMES gmp)
    find_library(GMPXX_LIBRARY NAMES gmpxx)
    if(GMP_INCLUDE_DIR AND GMP_LIBRARY AND GMPXX_LIBRARY)
        add_library(GMP::gmp UNKNOWN IMPORTED)
        set_target_properties(GMP::gmp PROPERTIES
            IMPORTED_LOCATION ${GMP_LIBRARY}
            INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR}
        )
        add_library(GMP::gmpxx UNKNOWN IMPORTED)
        set_target_properties(GMP::gmpxx PROPERTIES
            IMPORTED_LOCATION ${GMPXX_LIBRARY}
            INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR}
        )
    endif()
endif()

if(NOT TARGET MPFR::MPFR)
    find_path(MPFR_INCLUDE_DIR NAMES mpfr.h)
    find_library(MPFR_LIBRARY NAMES mpfr)
    if(MPFR_INCLUDE_DIR AND MPFR_LIBRARY)
        add_library(MPFR::MPFR UNKNOWN IMPORTED)
        set_target_properties(MPFR::MPFR PROPERTIES
            IMPORTED_LOCATION ${MPFR_LIBRARY}
            INTERFACE_INCLUDE_DIRECTORIES ${MPFR_INCLUDE_DIR}
        )
    endif()
endif()

# Determine availability of gridsynth dependencies and prompt if missing
set(NWQEC_MISSING_GRIDSYNTH_DEPS FALSE)
if(NOT (TARGET GMP::gmp AND TARGET GMP::gmpxx AND TARGET MPFR::MPFR))
    set(NWQEC_MISSING_GRIDSYNTH_DEPS TRUE)
endif()

if(NWQEC_MISSING_GRIDSYNTH_DEPS AND NOT NWQEC_ALLOW_NO_GMP)
    if(NOT Python3_Interpreter_FOUND)
        find_package(Python3 COMPONENTS Interpreter QUIET)
    endif()
    if(NOT Python3_Interpreter_FOUND)
        message(FATAL_ERROR "
*** GMP/MPFR MISSING ***\nInstall the GMP and MPFR development libraries, then rerun CMake.\n  macOS (Homebrew): brew install gmp mpfr\n  Ubuntu/Debian: sudo apt-get install -y libgmp-dev libmpfr-dev\n  Fedora: sudo dnf install gmp-devel mpfr-devel\n  Windows (MSYS2): pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr\n
If you prefer to proceed without the C++ gridsynth backend, rerun with\n  -DNWQEC_ALLOW_NO_GMP=ON\nor for pip:\n  pip install . --config-settings=cmake.define.NWQEC_ALLOW_NO_GMP=ON\n
")
    endif()

    set(_NWQEC_PROMPT_SCRIPT [=[
import sys
prompt = """GMP and MPFR were not detected.

Install instructions:
  macOS (Homebrew): brew install gmp mpfr
  Ubuntu/Debian: sudo apt-get install -y libgmp-dev libmpfr-dev
  Fedora: sudo dnf install gmp-devel mpfr-devel
  Windows (MSYS2): pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr

Install these libraries before continuing? [Y/n]: """
try:
    answer = input(prompt)
except EOFError:
    sys.exit(2)
answer = answer.strip().lower()
if not answer:
    answer = 'y'
print(answer)
]=])
    execute_process(COMMAND "${Python3_EXECUTABLE}" -c "${_NWQEC_PROMPT_SCRIPT}"
                    RESULT_VARIABLE _NWQEC_PROMPT_RESULT
                    OUTPUT_VARIABLE _NWQEC_PROMPT_OUTPUT)
    if(NOT _NWQEC_PROMPT_RESULT EQUAL 0)
        message(FATAL_ERROR "
*** GMP/MPFR MISSING ***\nConfiguration cancelled at your request. Install the libraries above or rerun with\n  -DNWQEC_ALLOW_NO_GMP=ON\n(pip: pip install . --config-settings=cmake.define.NWQEC_ALLOW_NO_GMP=ON)\n
")
    endif()

    string(STRIP "${_NWQEC_PROMPT_OUTPUT}" _NWQEC_PROMPT_OUTPUT)
    if(_NWQEC_PROMPT_OUTPUT STREQUAL "y" OR _NWQEC_PROMPT_OUTPUT STREQUAL "yes")
        message(FATAL_ERROR "Please install GMP and MPFR, then rerun configuration.")
    elseif(_NWQEC_PROMPT_OUTPUT STREQUAL "n" OR _NWQEC_PROMPT_OUTPUT STREQUAL "no")
        set(NWQEC_ALLOW_NO_GMP ON CACHE BOOL "Allow building without GMP/MPFR (use Python gridsynth fallback)" FORCE)
        message(WARNING "Proceeding without GMP/MPFR. Install python packages 'pygridsynth' and 'mpmath' to enable runtime RZ synthesis.")
    else()
        message(FATAL_ERROR "Unrecognized response: '${_NWQEC_PROMPT_OUTPUT}'. Aborting.")
    endif()
endif()

# Interface target for gridsynth that carries MP deps
add_library(nwqec_gridsynth INTERFACE)
target_link_libraries(nwqec_gridsynth INTERFACE nwqec)

set(_NWQEC_GRIDSYNTH_CPP_FLAG 0)
if(TARGET GMP::gmp AND TARGET GMP::gmpxx AND TARGET MPFR::MPFR)
    target_link_libraries(nwqec_gridsynth INTERFACE GMP::gmp GMP::gmpxx MPFR::MPFR)
    set(_NWQEC_GRIDSYNTH_CPP_FLAG 1)
endif()

target_compile_definitions(nwqec_gridsynth INTERFACE NWQEC_WITH_GRIDSYNTH_CPP=${_NWQEC_GRIDSYNTH_CPP_FLAG})

# =============================================================================
# Executables
# =============================================================================

add_executable(transpiler tools/transpiler.cpp)
set(_NWQEC_CLI_TARGETS transpiler)
target_link_libraries(transpiler PRIVATE nwqec_gridsynth)
target_compile_definitions(transpiler PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}" NWQEC_WITH_GRIDSYNTH_CPP=${_NWQEC_GRIDSYNTH_CPP_FLAG})
target_compile_options(transpiler PRIVATE ${COMMON_COMPILE_OPTIONS})

if(_NWQEC_GRIDSYNTH_CPP_FLAG)
    add_executable(gridsynth tools/gridsynth.cpp)
    target_link_libraries(gridsynth PRIVATE nwqec_gridsynth)
    target_compile_definitions(gridsynth PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}" NWQEC_WITH_GRIDSYNTH_CPP=${_NWQEC_GRIDSYNTH_CPP_FLAG})
    target_compile_options(gridsynth PRIVATE ${COMMON_COMPILE_OPTIONS})
    list(APPEND _NWQEC_CLI_TARGETS gridsynth)
endif()

# =============================================================================
# Build type specific tweaks
# =============================================================================

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        foreach(cli_target IN LISTS _NWQEC_CLI_TARGETS)
            if(NWQEC_ENABLE_NATIVE)
                target_compile_options(${cli_target} PRIVATE -O3 -march=native -funroll-loops)
            else()
                target_compile_options(${cli_target} PRIVATE -O3)
            endif()
        endforeach()
    else()
        foreach(cli_target IN LISTS _NWQEC_CLI_TARGETS)
            target_compile_options(${cli_target} PRIVATE -O2)
        endforeach()
    endif()

    if(NWQEC_ENABLE_LTO)
        include(CheckIPOSupported)
        check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_err)
        if(_ipo_ok)
            foreach(cli_target IN LISTS _NWQEC_CLI_TARGETS)
                set_property(TARGET ${cli_target} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
            endforeach()
        endif()
    endif()
endif()

message(STATUS "Configured executable: transpiler")
if(TARGET gridsynth)
    message(STATUS "Configured executable: gridsynth")
else()
    message(STATUS "Skipping gridsynth CLI (GMP/MPFR not available)")
endif()

# =============================================================================
# Install rules (CLI binaries and headers)
# =============================================================================
include(GNUInstallDirs)
install(TARGETS transpiler
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(TARGET gridsynth)
    install(TARGETS gridsynth RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# Install interface libraries and export targets
install(TARGETS nwqec nwqec_gridsynth
        EXPORT NWQECTargets)

install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT NWQECTargets
        NAMESPACE NWQEC::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NWQEC)

# Package config files
include(CMakePackageConfigHelpers)

# Record whether gridsynth deps were present at build
set(NWQEC_WITH_GRIDSYNTH OFF)
if(TARGET GMP::gmp AND TARGET GMP::gmpxx AND TARGET MPFR::MPFR)
    set(NWQEC_WITH_GRIDSYNTH ON)
endif()

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/NWQECConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/NWQECConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NWQEC
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/NWQECConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/NWQECConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/NWQECConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/NWQEC)

# =============================================================================
# Tests (CLI)
# =============================================================================
include(CTest)
if(BUILD_TESTING)
    add_test(NAME transpiler_help
             COMMAND $<TARGET_FILE:transpiler> --help)
    add_test(NAME transpiler_qft
             COMMAND $<TARGET_FILE:transpiler> --qft 4 --no-save)
    if(TARGET gridsynth)
        add_test(NAME gridsynth_basic
                 COMMAND $<TARGET_FILE:gridsynth> pi/4 10)
    endif()
endif()

# =============================================================================
# Python bindings (pybind11)
# =============================================================================
if(NWQEC_BUILD_PYTHON)
    find_package(pybind11 CONFIG QUIET)
    if(pybind11_FOUND)
        # Use a different CMake target name to avoid clashing with the INTERFACE library 'nwqec'
        pybind11_add_module(nwqec_ext MODULE python/nwqec/_core.cpp)
        # Ensure the produced module is packaged as nwqec._core
        set_target_properties(nwqec_ext PROPERTIES OUTPUT_NAME "_core")

        target_link_libraries(nwqec_ext PRIVATE nwqec)
        if(TARGET GMP::gmp AND TARGET GMP::gmpxx AND TARGET MPFR::MPFR)
            target_link_libraries(nwqec_ext PRIVATE GMP::gmp GMP::gmpxx MPFR::MPFR)
            target_compile_definitions(nwqec_ext PRIVATE NWQEC_WITH_GRIDSYNTH_CPP=1)
        else()
            target_compile_definitions(nwqec_ext PRIVATE NWQEC_WITH_GRIDSYNTH_CPP=0)
        endif()
        target_compile_definitions(nwqec_ext PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        target_compile_options(nwqec_ext PRIVATE ${COMMON_COMPILE_OPTIONS})
        # Install extension into the correct Python platlib dir when building a wheel
        if(DEFINED SKBUILD_PLATLIB_DIR)
            install(TARGETS nwqec_ext LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/nwqec)
        endif()
        message(STATUS "Configured Python module target: nwqec_ext (module name: nwqec._core)")
    else()
        message(WARNING "pybind11 not found; skipping Python module. Install with 'pip install pybind11' and reconfigure, or set NWQEC_BUILD_PYTHON=OFF.")
    endif()
endif()
