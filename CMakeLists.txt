cmake_minimum_required(VERSION 3.10)
project(ChemTrain VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directory for binary
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Sanitizer options
option(WITH_ASAN "Build with AddressSanitizer" OFF)
option(WITH_TSAN "Build with ThreadSanitizer" OFF)
option(WITH_UBSAN "Build with UndefinedBehaviorSanitizer" OFF)

# Configure sanitizers if enabled
if(WITH_ASAN)
    message(STATUS "Building with AddressSanitizer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

if(WITH_TSAN)
    message(STATUS "Building with ThreadSanitizer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=thread")
endif()

if(WITH_UBSAN)
    message(STATUS "Building with UndefinedBehaviorSanitizer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=undefined")
endif()

# RPATH configuration
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:/usr/local/lib:/usr/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Enable multithreading support
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)

# Define version information directly
add_compile_definitions(
    VERSION="${CMAKE_PROJECT_VERSION}"
    CHEMTRAIN_VERSION="${CMAKE_PROJECT_VERSION}"
)

# If find_package fails OR we skip it, locate libraries manually
message(STATUS "Forcing manual RDKit library search...")

# Set RDKit paths - try multiple common locations
set(RDKIT_INCLUDE_DIR "/usr/include/rdkit/")
set(RDKIT_LIBRARY_DIR "/usr/lib/x86_64-linux-gnu/")

# Add RDKit include directory
include_directories(${RDKIT_INCLUDE_DIR})

# Find RDKit libraries with more paths
set(RDKIT_LIB_SEARCH_PATHS
    ${RDKIT_LIBRARY_DIR}
    /usr/lib/x86_64-linux-gnu
    /usr/lib
    /usr/local/lib
    /lib
)

# Core libraries
find_library(RDKIT_GRAPHMOL_LIB NAMES GraphMol RDKitGraphMol PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_SMILEPARSE_LIB NAMES SmilesParse RDKitSmilesParse PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_CHEMTRANSFORMS_LIB NAMES ChemTransforms RDKitChemTransforms PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_DESCRIPTORS_LIB NAMES Descriptors RDKitDescriptors PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_RDGENERAL_LIB NAMES RDGeneral RDKitRDGeneral PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_MOLSTANDARDIZE_LIB NAMES MolStandardize RDKitMolStandardize PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_SUBSTRUCTMATCH_LIB NAMES SubstructMatch RDKitSubstructMatch PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_RDGEOMETRYLIB_LIB NAMES RDGeometryLib RDKitRDGeometryLib PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_DATASTRUCTS_LIB NAMES DataStructs RDKitDataStructs PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_FINGERPRINTS_LIB NAMES Fingerprints RDKitFingerprints PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_DEPICTOR_LIB NAMES Depictor RDKitDepictor PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_PARTIALCHARGES_LIB NAMES PartialCharges RDKitPartialCharges PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_MOLOPSLIB_LIB NAMES GraphMol RDKitGraphMol PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_RDLOGLIB_LIB NAMES RDGeneral RDKitRDGeneral PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_INVARLIB_LIB NAMES GraphMol RDKitGraphMol PATHS ${RDKIT_LIB_SEARCH_PATHS}) 

# Critical additions for your undefined references
find_library(RDKIT_MOLOPSLIB_LIB NAMES MolOps RDKitMolOps PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_RDLOGLIB_LIB NAMES RDLog RDKitRDLog PATHS ${RDKIT_LIB_SEARCH_PATHS})
find_library(RDKIT_INVARLIB_LIB NAMES Invar RDKitInvar PATHS ${RDKIT_LIB_SEARCH_PATHS})

# Set initial RDKIT_LIBRARIES list
set(RDKIT_LIBRARIES
    ${RDKIT_GRAPHMOL_LIB}
    ${RDKIT_SMILEPARSE_LIB}
    ${RDKIT_FILEPARSERS_LIB}
    ${RDKIT_CHEMTRANSFORMS_LIB}
    ${RDKIT_DESCRIPTORS_LIB}
    ${RDKIT_RDGENERAL_LIB}
    ${RDKIT_MOLSTANDARDIZE_LIB}
    ${RDKIT_SUBSTRUCTMATCH_LIB}
    ${RDKIT_RDGEOMETRYLIB_LIB}
    ${RDKIT_DATASTRUCTS_LIB}
    ${RDKIT_FINGERPRINTS_LIB}
    ${RDKIT_DEPICTOR_LIB}
    ${RDKIT_PARTIALCHARGES_LIB}
    ${RDKIT_MOLOPSLIB_LIB}
    ${RDKIT_RDLOGLIB_LIB}
    ${RDKIT_INVARLIB_LIB}
)

# Add the additional libraries if found
foreach(lib 
    RDKIT_BONDOPSLIB_LIB 
    RDKIT_ATOMOPSLIB_LIB
    RDKIT_SCAFFOLDNETWORK_LIB)
    if(${lib})
        list(APPEND RDKIT_LIBRARIES ${${lib}})
    endif()
endforeach()

# Additional commonly needed RDKit libraries
find_library(RDKIT_COORDGEN_LIB NAMES CoordGen RDKitCoordGen PATHS ${RDKIT_LIB_SEARCH_PATHS})
if(RDKIT_COORDGEN_LIB)
    list(APPEND RDKIT_LIBRARIES ${RDKIT_COORDGEN_LIB})
endif()

find_library(RDKIT_INCHI_LIB NAMES RDInchiLib PATHS ${RDKIT_LIB_SEARCH_PATHS})
if(RDKIT_INCHI_LIB)
    list(APPEND RDKIT_LIBRARIES ${RDKIT_INCHI_LIB})
endif()

# Find system libraries that RDKit might depend on
find_library(LOG_LIB NAMES log)
if(LOG_LIB)
    list(APPEND RDKIT_LIBRARIES ${LOG_LIB})
endif()

# Define RDKit manual flag directly
add_compile_definitions(WITH_RDKIT_MANUAL)

list(APPEND CMAKE_INSTALL_RPATH "${RDKIT_LIBRARY_DIR}")

# Add debug messages
message(STATUS "RDKIT_MOLOPSLIB_LIB: ${RDKIT_MOLOPSLIB_LIB}")
message(STATUS "RDKIT_SMILEPARSE_LIB: ${RDKIT_SMILEPARSE_LIB}")
message(STATUS "RDKIT_GRAPHMOL_LIB: ${RDKIT_GRAPHMOL_LIB}")

# Find Boost (required by RDKit)
find_package(Boost REQUIRED COMPONENTS log log_setup thread filesystem system date_time regex atomic)
include_directories(${Boost_INCLUDE_DIRS})

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/src
)

# Source files
set(SOURCES
src/common.cpp
    src/registry.cpp
    src/descriptors/regex_utils.cpp
    src/descriptors/regexfractions.cpp
    src/descriptors/weightfractions.cpp
    src/descriptors/regexelectronegativity.cpp
    src/descriptors/bordwellscale.cpp
    src/descriptors/bonds.cpp
    src/descriptors/radius.cpp
    src/descriptors/brotomoreau.cpp
    src/descriptors/graphoperations.cpp
    src/descriptors/c_wrappers.cpp
)

# C source files
set(C_SOURCES
    src/csv_parser.c
    src/csv_writer.c
    src/descriptors/topological.c
    src/descriptors/bitsoperations.c
    src/descriptors/regexcounts.c
)

# Add CSV parser support
add_compile_definitions(WITH_FAST_CSV_PARSER)

# Add script to run before build
add_custom_target(update_registry
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/update_registry.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Updating descriptor registry"
)

# Set C standard for C files
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Create object library for C files
add_library(csv_utils OBJECT ${C_SOURCES})
set_target_properties(csv_utils PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_STANDARD 11
    C_STANDARD_REQUIRED ON
)

# Executable
add_executable(chemtrain ${SOURCES} $<TARGET_OBJECTS:csv_utils> src/main.cpp)

# Make the main target depend on the registry update
add_dependencies(chemtrain update_registry)

# Include directories
target_include_directories(chemtrain PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Find RE2 library
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(RE2 IMPORTED_TARGET re2)
endif()

if(NOT RE2_FOUND)
    find_library(RE2_LIBRARY NAMES re2)
    find_path(RE2_INCLUDE_DIR NAMES re2/re2.h)
    
    if(RE2_LIBRARY AND RE2_INCLUDE_DIR)
        set(RE2_FOUND TRUE)
        set(RE2_LIBRARIES ${RE2_LIBRARY})
        set(RE2_INCLUDE_DIRS ${RE2_INCLUDE_DIR})
        message(STATUS "Found RE2 manually: ${RE2_LIBRARY}")
    else()
        message(FATAL_ERROR "RE2 library not found. Install libre2-dev package or equivalent")
    endif()
endif()

# Include RE2 in the build
include_directories(${RE2_INCLUDE_DIRS})

# Replace the existing Arrow/Parquet detection
find_package(PkgConfig REQUIRED)
pkg_check_modules(ARROW REQUIRED arrow>=20.0.0)
pkg_check_modules(PARQUET REQUIRED parquet>=20.0.0)

if(ARROW_FOUND AND PARQUET_FOUND)
    add_compile_definitions(WITH_ARROW_PARQUET)
    message(STATUS "Found Arrow and Parquet, enabling Parquet file support")
    include_directories(${ARROW_INCLUDE_DIRS} ${PARQUET_INCLUDE_DIRS})
endif()

# Add OpenMP support
find_package(OpenMP REQUIRED)

# Final target_link_libraries call
target_link_libraries(chemtrain PRIVATE
    ${RDKIT_LIBRARIES}
    ${Boost_LIBRARIES}
    ${Boost_LOG_LIBRARY}
    ${Boost_LOG_SETUP_LIBRARY}
    ${Boost_THREAD_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_DATE_TIME_LIBRARY}
    ${Boost_REGEX_LIBRARY}
    ${Boost_ATOMIC_LIBRARY}
    ${RE2_LIBRARIES}
    Threads::Threads
    ${ARROW_LIBRARIES}
    ${PARQUET_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(chemtrain PRIVATE OpenMP::OpenMP_CXX)
endif()

# Set output name
set_target_properties(chemtrain PROPERTIES
    OUTPUT_NAME "chemtrain"
    BUILD_WITH_INSTALL_RPATH FALSE
    SKIP_BUILD_RPATH FALSE
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:/usr/local/lib:/usr/lib:${RDKIT_LIBRARY_DIR}"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

# Create a config.hpp file directly in the build directory
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/src/config.hpp 
"#pragma once

// Version information
#define VERSION \"${CMAKE_PROJECT_VERSION}\"
#define CHEMTRAIN_VERSION \"${CMAKE_PROJECT_VERSION}\"

// Feature flags - automatically defined based on find results
#ifdef WITH_FAST_CSV_PARSER
// Fast CSV parser is available
#endif

#ifdef WITH_RDKIT_MANUAL
// Using manual RDKit integration
#endif

#ifdef WITH_ARROW_PARQUET
// Arrow/Parquet support is available
#endif

// Sanitizer flags
#ifdef WITH_ASAN
// Address Sanitizer is enabled
#endif

#ifdef WITH_TSAN
// Thread Sanitizer is enabled
#endif

#ifdef WITH_UBSAN
// Undefined Behavior Sanitizer is enabled
#endif
")

include_directories(${CMAKE_CURRENT_BINARY_DIR}/src)

# Install targets
install(TARGETS chemtrain
        RUNTIME DESTINATION bin)

# Install headers
install(FILES
    src/common.hpp
    src/registry.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/src/config.hpp
    DESTINATION include/chemtrain)

# After build, copy chemtrain binary to bin/ directory
add_custom_command(TARGET chemtrain POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:chemtrain> ${CMAKE_SOURCE_DIR}/bin/
    COMMENT "Copying chemtrain binary to bin/ directory"
)