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

# Directory for downloading prebuilt GMP/MPFR binaries
set(NWQEC_PREBUILT_DIR "${CMAKE_BINARY_DIR}/prebuilt_deps" CACHE PATH "Download cache for prebuilt dependencies")

# Get prebuilt binary info for current platform
function(nwqec_get_prebuilt_info out_url out_hash out_name)
    set(_url "")
    set(_hash "")
    set(_name "")

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64)$")
        set(_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-linux-x86_64.zip")
        set(_hash "b809417ec3ac3720bf8f359a1b7385197e830a2687a92ac60f8c5d013f412633")
        set(_name "gmp-mpfr-linux-x86_64.zip")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
        set(_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-linux-arm64.zip")
        set(_hash "48b1e9f01659c9072445b6430f2488467d88fc456cddfb7e423787a4168f8bf4")
        set(_name "gmp-mpfr-linux-arm64.zip")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
        set(_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-macos-arm64.zip")
        set(_hash "b247e76817ea3b0e7ed98b918c893ff8db301db48fcf52bba9f84412312b1f53")
        set(_name "gmp-mpfr-macos-arm64.zip")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64)$")
        set(_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-macos-x86_64.zip")
        set(_hash "5fcc8af345ade29dabb48e86a94852c66f3926d6bf82209cdd94e73725dcd1be")
        set(_name "gmp-mpfr-macos-x86_64.zip")
    endif()

    set(${out_url} "${_url}" PARENT_SCOPE)
    set(${out_hash} "${_hash}" PARENT_SCOPE)
    set(${out_name} "${_name}" PARENT_SCOPE)
endfunction()

function(nwqec_pick_library_from_bundle out_var lib_dir base_name)
    set(_candidates)
    if(CMAKE_STATIC_LIBRARY_PREFIX)
        list(APPEND _candidates "${lib_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}${base_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
    if(CMAKE_SHARED_LIBRARY_PREFIX)
        list(APPEND _candidates "${lib_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}${base_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
    list(APPEND _candidates
        "${lib_dir}/lib${base_name}.a"
        "${lib_dir}/lib${base_name}.lib"
        "${lib_dir}/${base_name}.a"
        "${lib_dir}/${base_name}.lib")
    if(APPLE)
        list(APPEND _candidates "${lib_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}${base_name}.dylib")
    endif()
    if(WIN32)
        list(APPEND _candidates "${lib_dir}/${base_name}.lib" "${lib_dir}/${base_name}.dll")
    endif()
    list(REMOVE_DUPLICATES _candidates)

    foreach(_candidate IN LISTS _candidates)
        if(_candidate AND EXISTS "${_candidate}")
            set(${out_var} "${_candidate}" PARENT_SCOPE)
            return()
        endif()
    endforeach()

    file(GLOB _versioned_candidates
        "${lib_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}${base_name}.*"
        "${lib_dir}/lib${base_name}.*"
        "${lib_dir}/${base_name}.*")
    if(_versioned_candidates)
        list(FILTER _versioned_candidates INCLUDE REGEX "\\.(a|so(\\.[0-9]+)*|dylib|dll|lib)$")
        list(SORT _versioned_candidates)
        list(GET _versioned_candidates 0 _first_candidate)
        set(${out_var} "${_first_candidate}" PARENT_SCOPE)
        return()
    endif()

    set(${out_var} "" PARENT_SCOPE)
endfunction()

# Download and setup prebuilt GMP/MPFR binaries
function(nwqec_setup_prebuilt_deps out_root)
    nwqec_get_prebuilt_info(_url _hash _filename)
    if(_url STREQUAL "")
        set(${out_root} "" PARENT_SCOPE)
        return()
    endif()

    file(MAKE_DIRECTORY "${NWQEC_PREBUILT_DIR}")
    set(_zip_path "${NWQEC_PREBUILT_DIR}/${_filename}")
    
    # Check if we need to download
    set(_need_download TRUE)
    if(EXISTS "${_zip_path}")
        file(SHA256 "${_zip_path}" _existing_hash)
        if(_existing_hash STREQUAL "${_hash}")
            set(_need_download FALSE)
        else()
            message(STATUS "Cached file hash mismatch; re-downloading ${_filename}")
            file(REMOVE "${_zip_path}")
        endif()
    endif()

    if(_need_download)
        message(STATUS "Downloading prebuilt GMP/MPFR from ${_url}")
        file(DOWNLOAD "${_url}" "${_zip_path}" SHOW_PROGRESS EXPECTED_HASH SHA256=${_hash})
    endif()

    # Extract to platform-specific directory
    string(REGEX REPLACE "\\.zip$" "" _extract_dir "${_filename}")
    set(_extract_path "${NWQEC_PREBUILT_DIR}/${_extract_dir}")
    set(_marker_file "${_extract_path}/.extracted")
    
    if(NOT EXISTS "${_marker_file}")
        if(EXISTS "${_extract_path}")
            file(REMOVE_RECURSE "${_extract_path}")
        endif()
        file(MAKE_DIRECTORY "${_extract_path}")
        
        # Extract zip -> find tarball -> extract tarball
        execute_process(
            COMMAND "${CMAKE_COMMAND}" -E tar xf "${_zip_path}"
            WORKING_DIRECTORY "${_extract_path}"
            RESULT_VARIABLE _zip_result)
        if(NOT _zip_result EQUAL 0)
            message(FATAL_ERROR "Failed to extract ${_zip_path}")
        endif()

        # Find and extract the tarball inside
        file(GLOB _tarballs "${_extract_path}/*.tar.gz")
        if(NOT _tarballs)
            message(FATAL_ERROR "No .tar.gz found in ${_zip_path}")
        endif()
        list(GET _tarballs 0 _tarball)

        execute_process(
            COMMAND "${CMAKE_COMMAND}" -E tar xzf "${_tarball}"
            WORKING_DIRECTORY "${_extract_path}"
            RESULT_VARIABLE _tar_result)
        if(NOT _tar_result EQUAL 0)
            message(FATAL_ERROR "Failed to extract ${_tarball}")
        endif()

        file(WRITE "${_marker_file}" "extracted")
    endif()

    set(${out_root} "${_extract_path}" PARENT_SCOPE)
endfunction()

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
if(WIN32)
    target_compile_definitions(nwqec INTERFACE _USE_MATH_DEFINES)
endif()

# =============================================================================
# Gridsynth dependency (GMP/MPFR) and target
# =============================================================================

if(WIN32)
    # Windows: Skip GMP/MPFR due to complexity, use pygridsynth fallback
    message(STATUS "Windows detected: skipping C++ gridsynth (use pygridsynth fallback)")
    set(NWQEC_WITH_GRIDSYNTH OFF)
else()
    # macOS/Linux: Always use prebuilt GMP/MPFR binaries
    message(STATUS "Setting up prebuilt GMP/MPFR for ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
    nwqec_setup_prebuilt_deps(_PREBUILT_ROOT)

    if(NOT _PREBUILT_ROOT)
        message(FATAL_ERROR "No prebuilt GMP/MPFR available for ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
    endif()

    set(_INCLUDE_DIR "${_PREBUILT_ROOT}/include")
    set(_LIB_DIR "${_PREBUILT_ROOT}/lib")
    
    if(NOT EXISTS "${_INCLUDE_DIR}" OR NOT EXISTS "${_LIB_DIR}")
        message(FATAL_ERROR "Prebuilt archive missing include/lib directories at ${_PREBUILT_ROOT}")
    endif()

    # Find the library files
    nwqec_pick_library_from_bundle(_GMP_LIB "${_LIB_DIR}" gmp)
    nwqec_pick_library_from_bundle(_MPFR_LIB "${_LIB_DIR}" mpfr)

    if(NOT _GMP_LIB OR NOT _MPFR_LIB)
        message(FATAL_ERROR "Required libraries not found in ${_LIB_DIR}")
    endif()

    # Create imported targets
    add_library(GMP::gmp UNKNOWN IMPORTED)
    set_target_properties(GMP::gmp PROPERTIES
        IMPORTED_LOCATION "${_GMP_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${_INCLUDE_DIR}"
    )

    add_library(MPFR::MPFR UNKNOWN IMPORTED)
    set_target_properties(MPFR::MPFR PROPERTIES
        IMPORTED_LOCATION "${_MPFR_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${_INCLUDE_DIR}"
    )

    set(NWQEC_WITH_GRIDSYNTH ON)
    message(STATUS "Using prebuilt GMP/MPFR from ${_PREBUILT_ROOT}")
endif()

# Interface target for gridsynth that carries MP deps (only if available)
if(NWQEC_WITH_GRIDSYNTH)
    add_library(nwqec_gridsynth INTERFACE)
    target_link_libraries(nwqec_gridsynth INTERFACE nwqec GMP::gmp MPFR::MPFR)
    target_compile_definitions(nwqec_gridsynth INTERFACE NWQEC_WITH_GRIDSYNTH_CPP=1)
else()
    add_library(nwqec_gridsynth INTERFACE)
    target_link_libraries(nwqec_gridsynth INTERFACE nwqec)
    target_compile_definitions(nwqec_gridsynth INTERFACE NWQEC_WITH_GRIDSYNTH_CPP=0)
endif()

# =============================================================================
# Executables
# =============================================================================

if(NWQEC_WITH_GRIDSYNTH)
    add_executable(nwqec-cli tools/nwqec.cpp)
    set(_NWQEC_CLI_TARGETS nwqec-cli)
    target_link_libraries(nwqec-cli PRIVATE nwqec_gridsynth)
    target_compile_definitions(nwqec-cli PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    target_compile_options(nwqec-cli PRIVATE ${COMMON_COMPILE_OPTIONS})

    add_executable(gridsynth tools/gridsynth.cpp)
    target_link_libraries(gridsynth PRIVATE nwqec_gridsynth)
    target_compile_definitions(gridsynth PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    target_compile_options(gridsynth PRIVATE ${COMMON_COMPILE_OPTIONS})
    list(APPEND _NWQEC_CLI_TARGETS gridsynth)
else()
    message(STATUS "Skipping CLI tools (GMP/MPFR not available on this platform)")
    message(STATUS "  → Use Python bindings instead: import nwqec")
    set(_NWQEC_CLI_TARGETS)
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

# Status messages
if(NWQEC_WITH_GRIDSYNTH)
    message(STATUS "Configured executable: nwqec-cli (with GMP/MPFR gridsynth backend)")
    message(STATUS "Configured executable: gridsynth (with GMP/MPFR)")
else()
    message(STATUS "Skipped CLI tools (C++ gridsynth unavailable)")
    message(STATUS "  → Use Python bindings instead: import nwqec")
    message(STATUS "  → For RZ synthesis, install pygridsynth: pip install pygridsynth")
endif()

# =============================================================================
# Install rules (CLI binaries and headers)
# =============================================================================
include(GNUInstallDirs)
if(NWQEC_WITH_GRIDSYNTH)
    install(TARGETS nwqec-cli gridsynth
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
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

# Record whether gridsynth deps were present at build (already set above)

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
if(BUILD_TESTING AND NWQEC_WITH_GRIDSYNTH)
    add_test(NAME nwqec_help
             COMMAND $<TARGET_FILE:nwqec-cli> --help)
    add_test(NAME nwqec_qft
             COMMAND $<TARGET_FILE:nwqec-cli> --qft 4 --no-save)
    add_test(NAME gridsynth_basic
             COMMAND $<TARGET_FILE:gridsynth> pi/4 10)
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

        target_link_libraries(nwqec_ext PRIVATE nwqec_gridsynth)
        target_compile_definitions(nwqec_ext PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        
        # Explicitly set gridsynth availability (redundant but ensures it's set)
        if(NWQEC_WITH_GRIDSYNTH)
            target_compile_definitions(nwqec_ext PRIVATE NWQEC_WITH_GRIDSYNTH_CPP=1)
        else()
            target_compile_definitions(nwqec_ext PRIVATE NWQEC_WITH_GRIDSYNTH_CPP=0)
        endif()
        
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
