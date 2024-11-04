# cmake-format: off
#
# Variables to control this CMake build:
# ===================================================================================
#    CMake flag       |  environment variable  |  Description
# ===================================================================================
# -DCMAKE_CXX_STANDARD|                        | standard for C++ compilation
# ===================================================================================
#
#
# ==============================================================
#    DEPENDENCIES (used by find_package)
# ==============================================================
#      Name           |  Description
# ==============================================================
# EDM4HEP             |   Event data model based on podio
# podio               |   IO library
# DD4hep              |   Geometry framework
# ROOT                |   CERN ROOT
# EDM4EIC
# ==============================================================
#
# cmake-format: on

cmake_minimum_required(VERSION 3.24)

project(recoEvents
	DESCRIPTION "Project the events TTree of the podio_output of eicrecon")

# Set version based on git (not working as of 07/03: fatal: No names found, cannot describe anything.)
#include(cmake/git_version.cmake)
#set_git_version(CMAKE_PROJECT_VERSION)

# Set default standard to C++20
set(CMAKE_CXX_STANDARD_MIN 20)
set(CMAKE_CXX_STANDARD
    ${CMAKE_CXX_STANDARD_MIN}
    CACHE STRING "C++ standard to be used")
if(CMAKE_CXX_STANDARD LESS CMAKE_CXX_STANDARD_MIN)
  message(
    FATAL_ERROR
      "Unsupported C++ standard: ${CMAKE_CXX_STANDARD} (at least ${CMAKE_CXX_STANDARD_MIN} required)"
  )
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Debug)

### where to find the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

### add CERN ROOT
find_package(ROOT 6 REQUIRED)
include(${ROOT_USE_FILE})

# PODIO, EDM4HEP, EDM4EIC event models
find_package(podio REQUIRED)
find_package(EDM4HEP 0.7.1 REQUIRED)
message(STATUS "EDM4HEP_ROOT ${EDM4HEP_FOUND} ${EDM4HEP_INCLUDE_DIR}")
find_package(EDM4EIC 5.0 REQUIRED)
message(STATUS "EDM4EIC_ROOT ${EDM4EIC_FOUND} ${EDM4EIC_INCLUDE_DIR}")

find_package(DD4hep REQUIRED)
message(STATUS "DD4HEP_ROOT ${DD4hep_FOUND} ${DD4hep_INCLUDE_DIRS} ${DD4hep_LIB}")

include_directories(${CMAKE_SOURCE_DIR}/src /opt/local/include)
add_subdirectory(src)

# Install to the top directory by default
if(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
  set(CMAKE_INSTALL_PREFIX
      ${CMAKE_SOURCE_DIR}
      CACHE PATH "Install in top directory by default" FORCE)
endif()

# Install to standard location
include(GNUInstallDirs)

install(TARGETS recoEvents DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES "${PROJECT_BINARY_DIR}/src/librecoEventsDict_rdict.pcm" DESTINATION ${CMAKE_INSTALL_PREFIX})
