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
option(NWQEC_FORCE_BUNDLED_GMP_MPFR "Always download and use the bundled GMP/MPFR distribution" OFF)

set(NWQEC_GMP_MPFR_BUNDLE_URL_OVERRIDE "" CACHE STRING "Custom URL for a GMP/MPFR prebuilt bundle (.zip)")
set(NWQEC_GMP_MPFR_BUNDLE_HASH_OVERRIDE "" CACHE STRING "SHA256 for the overridden GMP/MPFR bundle")
set(NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/nwqec_gmp_mpfr" CACHE PATH "Download cache location for GMP/MPFR bundles")

mark_as_advanced(
    NWQEC_GMP_MPFR_BUNDLE_URL_OVERRIDE
    NWQEC_GMP_MPFR_BUNDLE_HASH_OVERRIDE
    NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR)

function(nwqec_select_gmp_mpfr_bundle out_url out_hash out_label)
    if(NOT "${NWQEC_GMP_MPFR_BUNDLE_URL_OVERRIDE}" STREQUAL "")
        set(${out_url} "${NWQEC_GMP_MPFR_BUNDLE_URL_OVERRIDE}" PARENT_SCOPE)
        set(${out_hash} "${NWQEC_GMP_MPFR_BUNDLE_HASH_OVERRIDE}" PARENT_SCOPE)
        if("${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "")
            set(_override_label "override")
        else()
            set(_override_label "override-${CMAKE_HOST_SYSTEM_NAME}")
        endif()
        set(${out_label} "${_override_label}" PARENT_SCOPE)
        return()
    endif()

    set(_bundle_url "")
    set(_bundle_hash "")
    set(_bundle_label "")

    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64)$")
        set(_bundle_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-linux-x86_64.zip")
        set(_bundle_hash "e4666dc44e3907da82c39a3cb91a4206cb24d43af4334bb03f56740f8d5a166b")
        set(_bundle_label "linux-x86_64")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
        set(_bundle_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-linux-arm64.zip")
        set(_bundle_hash "c8175833f8fa67fa521fc9b397f0f02c74cc39e7bdc124abfaae6f4a96dc16d5")
        set(_bundle_label "linux-arm64")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
        set(_bundle_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-macos-arm64.zip")
        set(_bundle_hash "775b6c7232ea46b0d245952aacd1616fed6510c9330dcf753a6db70f8f083fea")
        set(_bundle_label "macos-arm64")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64)$")
        set(_bundle_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-macos-x86_64.zip")
        set(_bundle_hash "10b253b7f075189d09c1fd52296977f31a45d27ac6bd5fe64a65d7ef4ac3107f")
        set(_bundle_label "macos-x86_64")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows" AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
        set(_bundle_url "https://github.com/pnnl/nwqec/releases/download/v0.1.0/gmp-mpfr-windows-x86_64.zip")
        set(_bundle_hash "7f407b22df082dfa942329e09b3bc25cdebfed3577486285352767ab5245586a")
        set(_bundle_label "windows-x86_64")
    endif()

    set(${out_url} "${_bundle_url}" PARENT_SCOPE)
    set(${out_hash} "${_bundle_hash}" PARENT_SCOPE)
    set(${out_label} "${_bundle_label}" PARENT_SCOPE)
endfunction()

function(nwqec_pick_library_from_bundle out_var lib_dir base_name)
    set(_candidates)
    if(CMAKE_STATIC_LIBRARY_PREFIX)
        list(APPEND _candidates "${lib_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}${base_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
    if(CMAKE_SHARED_LIBRARY_PREFIX)
        list(APPEND _candidates "${lib_dir}/${CMAKE_SHARED_LIBRARY_PREFIX}${base_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
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

function(nwqec_prepare_bundled_gmp_mpfr out_root)
    nwqec_select_gmp_mpfr_bundle(_bundle_url _bundle_hash _bundle_label)
    if(_bundle_url STREQUAL "")
        set(${out_root} "" PARENT_SCOPE)
        return()
    endif()

    file(MAKE_DIRECTORY "${NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR}")

    get_filename_component(_zip_filename "${_bundle_url}" NAME)
    if(_zip_filename STREQUAL "")
        set(_zip_filename "gmp-mpfr-bundle.zip")
    endif()
    set(_zip_path "${NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR}/${_zip_filename}")

    if(EXISTS "${_zip_path}" AND NOT _bundle_hash STREQUAL "")
        file(SHA256 "${_zip_path}" _existing_hash)
        if(NOT _existing_hash STREQUAL "${_bundle_hash}")
            message(STATUS "Cached GMP/MPFR bundle has mismatched hash; re-downloading ${_zip_filename}")
            file(REMOVE "${_zip_path}")
        endif()
    endif()

    if(NOT EXISTS "${_zip_path}")
        message(STATUS "Downloading GMP/MPFR bundle (${_bundle_label}) from ${_bundle_url}")
        if(_bundle_hash STREQUAL "")
            file(DOWNLOAD "${_bundle_url}" "${_zip_path}" SHOW_PROGRESS)
        else()
            file(DOWNLOAD "${_bundle_url}" "${_zip_path}" SHOW_PROGRESS EXPECTED_HASH SHA256=${_bundle_hash})
        endif()
    endif()

    string(REGEX REPLACE "\\.zip$" "" _bundle_base "${_zip_filename}")
    if(_bundle_base STREQUAL "${_zip_filename}")
        # If the archive name does not end in .zip (unexpected), add a suffix to avoid collisions.
        set(_bundle_base "${_bundle_base}-unzipped")
    endif()

    set(_tarball_path "")
    file(GLOB _existing_tarballs "${NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR}/*.tar.gz")
    foreach(_candidate IN LISTS _existing_tarballs)
        if(_bundle_label AND _candidate MATCHES "${_bundle_label}")
            set(_tarball_path "${_candidate}")
            break()
        elseif(NOT _bundle_label AND _candidate MATCHES "${_bundle_base}")
            set(_tarball_path "${_candidate}")
            break()
        endif()
    endforeach()
    if(_tarball_path STREQUAL "" AND _existing_tarballs)
        list(GET _existing_tarballs 0 _tarball_path)
    endif()

    if(_tarball_path STREQUAL "")
        execute_process(
            COMMAND "${CMAKE_COMMAND}" -E tar xvf "${_zip_path}"
            WORKING_DIRECTORY "${NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR}"
            RESULT_VARIABLE _zip_result
            OUTPUT_VARIABLE _zip_output
            ERROR_VARIABLE _zip_error)
        if(NOT _zip_result EQUAL 0)
            message(FATAL_ERROR "Failed to extract ${_zip_path} (exit code ${_zip_result}): ${_zip_error}")
        endif()

        file(GLOB _extracted_tarballs "${NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR}/*.tar.gz")
        foreach(_candidate IN LISTS _extracted_tarballs)
            if(_bundle_label AND _candidate MATCHES "${_bundle_label}")
                set(_tarball_path "${_candidate}")
                break()
            elseif(NOT _bundle_label AND _candidate MATCHES "${_bundle_base}")
                set(_tarball_path "${_candidate}")
                break()
            endif()
        endforeach()
        if(_tarball_path STREQUAL "" AND _extracted_tarballs)
            list(GET _extracted_tarballs 0 _tarball_path)
        endif()
    endif()

    if(_tarball_path STREQUAL "")
        message(FATAL_ERROR "Bundle ${_zip_filename} did not contain a .tar.gz payload that matches ${_bundle_label} or ${_bundle_base}.")
    endif()

    set(_bundle_prefix "${NWQEC_GMP_MPFR_BUNDLE_DOWNLOAD_DIR}/${_bundle_base}")
    set(_bundle_marker "${_bundle_prefix}/.nwqec-ready")
    if(NOT EXISTS "${_bundle_marker}")
        if(EXISTS "${_bundle_prefix}")
            file(REMOVE_RECURSE "${_bundle_prefix}")
        endif()
        file(MAKE_DIRECTORY "${_bundle_prefix}")
        execute_process(
            COMMAND "${CMAKE_COMMAND}" -E tar xzf "${_tarball_path}"
            WORKING_DIRECTORY "${_bundle_prefix}"
            RESULT_VARIABLE _tar_result)
        if(NOT _tar_result EQUAL 0)
            message(FATAL_ERROR "Failed to extract ${_tarball_path} (exit code ${_tar_result})")
        endif()
        file(WRITE "${_bundle_marker}" "ready")
    endif()

    set(${out_root} "${_bundle_prefix}" PARENT_SCOPE)
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

# Try find_package first; if unavailable, create imported targets from paths
if(NOT NWQEC_FORCE_BUNDLED_GMP_MPFR)
    find_package(GMP QUIET)
    find_package(MPFR QUIET)

    if(NOT TARGET GMP::gmp)
        find_path(GMP_INCLUDE_DIR NAMES gmp.h)
        find_library(GMP_LIBRARY NAMES gmp)
        if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
            add_library(GMP::gmp UNKNOWN IMPORTED)
            set_target_properties(GMP::gmp PROPERTIES
                IMPORTED_LOCATION ${GMP_LIBRARY}
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
endif()

if(NOT (TARGET GMP::gmp AND TARGET MPFR::MPFR))
    nwqec_prepare_bundled_gmp_mpfr(_NWQEC_GMP_MPFR_ROOT)

    if(_NWQEC_GMP_MPFR_ROOT)
        set(_NWQEC_GMP_MPFR_INCLUDE_DIR "${_NWQEC_GMP_MPFR_ROOT}/include")
        set(_NWQEC_GMP_MPFR_LIB_DIR "${_NWQEC_GMP_MPFR_ROOT}/lib")
        if(NOT EXISTS "${_NWQEC_GMP_MPFR_INCLUDE_DIR}" OR NOT EXISTS "${_NWQEC_GMP_MPFR_LIB_DIR}")
            message(FATAL_ERROR "Bundled GMP/MPFR archive is missing expected include/lib directories at ${_NWQEC_GMP_MPFR_ROOT}")
        endif()

        nwqec_pick_library_from_bundle(_NWQEC_GMP_LIB "${_NWQEC_GMP_MPFR_LIB_DIR}" gmp)
        nwqec_pick_library_from_bundle(_NWQEC_MPFR_LIB "${_NWQEC_GMP_MPFR_LIB_DIR}" mpfr)

        if(NOT _NWQEC_GMP_LIB OR NOT _NWQEC_MPFR_LIB)
            message(FATAL_ERROR "Bundled GMP/MPFR archive at ${_NWQEC_GMP_MPFR_ROOT} is missing required libraries.")
        endif()

        if(NOT TARGET GMP::gmp)
            add_library(GMP::gmp UNKNOWN IMPORTED)
            set_target_properties(GMP::gmp PROPERTIES
                IMPORTED_LOCATION "${_NWQEC_GMP_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${_NWQEC_GMP_MPFR_INCLUDE_DIR}"
            )
        endif()

        if(NOT TARGET MPFR::MPFR)
            add_library(MPFR::MPFR UNKNOWN IMPORTED)
            set_target_properties(MPFR::MPFR PROPERTIES
                IMPORTED_LOCATION "${_NWQEC_MPFR_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${_NWQEC_GMP_MPFR_INCLUDE_DIR}"
            )
        endif()

        set(_NWQEC_USING_BUNDLED_GMP_MPFR TRUE)
        message(STATUS "Using bundled GMP/MPFR package from ${_NWQEC_GMP_MPFR_ROOT}")
    endif()
endif()

if(NOT (TARGET GMP::gmp AND TARGET MPFR::MPFR))
    message(FATAL_ERROR "GMP/MPFR libraries are required but were not found, and no compatible bundled archive is defined for this platform. Set NWQEC_GMP_MPFR_BUNDLE_URL_OVERRIDE (and hash) or install GMP/MPFR development packages manually.")
endif()

# Interface target for gridsynth that carries MP deps
add_library(nwqec_gridsynth INTERFACE)
target_link_libraries(nwqec_gridsynth INTERFACE nwqec GMP::gmp MPFR::MPFR)

# =============================================================================
# Executables
# =============================================================================

add_executable(transpiler tools/transpiler.cpp)
set(_NWQEC_CLI_TARGETS transpiler)
target_link_libraries(transpiler PRIVATE nwqec_gridsynth)
target_compile_definitions(transpiler PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(transpiler PRIVATE ${COMMON_COMPILE_OPTIONS})

add_executable(gridsynth tools/gridsynth.cpp)
target_link_libraries(gridsynth PRIVATE nwqec_gridsynth)
target_compile_definitions(gridsynth PRIVATE PROJECT_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_options(gridsynth PRIVATE ${COMMON_COMPILE_OPTIONS})
list(APPEND _NWQEC_CLI_TARGETS gridsynth)

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
message(STATUS "Configured executable: gridsynth")

# =============================================================================
# Install rules (CLI binaries and headers)
# =============================================================================
include(GNUInstallDirs)
install(TARGETS transpiler gridsynth
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

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
set(NWQEC_WITH_GRIDSYNTH ON)

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

        target_link_libraries(nwqec_ext PRIVATE nwqec GMP::gmp MPFR::MPFR)
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
