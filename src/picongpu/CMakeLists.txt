# Copyright 2013-2017 Axel Huebl, Benjamin Schneider, Felix Schmitt, Heiko Burau,
#                     Rene Widera, Alexander Grund, Alexander Matthes
#
# This file is part of PIConGPU.
#
# PIConGPU is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIConGPU is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PIConGPU.
# If not, see <http://www.gnu.org/licenses/>.
#

################################################################################
# Required cmake version
################################################################################

cmake_minimum_required(VERSION 3.3.0)


################################################################################
# Project
################################################################################

project(PIConGPU)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "install prefix" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# set helper pathes to find libraries and packages
# Add specific hints
list(APPEND CMAKE_PREFIX_PATH "$ENV{MPI_ROOT}")
list(APPEND CMAKE_PREFIX_PATH "$ENV{CUDA_ROOT}")
list(APPEND CMAKE_PREFIX_PATH "$ENV{BOOST_ROOT}")
list(APPEND CMAKE_PREFIX_PATH "$ENV{HDF5_ROOT}")
list(APPEND CMAKE_PREFIX_PATH "$ENV{ADIOS_ROOT}")
# Add from environment after specific env vars
list(APPEND CMAKE_PREFIX_PATH "$ENV{CMAKE_PREFIX_PATH}")


###############################################################################
# Language Flags
###############################################################################

# enforce C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 11)


################################################################################
# Disallow in-source build
################################################################################

get_filename_component(SOURCE_DIR_ROOT ${PIConGPU_SOURCE_DIR}/../.. ABSOLUTE)
string(FIND "${PIConGPU_BINARY_DIR}"
            "${SOURCE_DIR_ROOT}" IN_SRC_POS)
if(IN_SRC_POS GREATER -1)
  message(FATAL_ERROR
    "PICoNGPU requires an out of source build. "
    "Please remove \n"
    "  - CMakeCache.txt\n"
    "  - CMakeFiles/\n"
    "and create a separate build directory. "
    "See: INSTALL.rst")
endif()

unset(IN_SRC_POS)


################################################################################
# Find PMacc
################################################################################

find_package(PMacc REQUIRED CONFIG PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../libPMacc")
include_directories(SYSTEM ${PMacc_INCLUDE_DIRS})
set(LIBS ${LIBS} ${PMacc_LIBRARIES})
add_definitions(${PMacc_DEFINITIONS})


################################################################################
# Find CUDA
################################################################################
find_package(CUDA 7.5 REQUIRED)


################################################################################
# Find MPI
################################################################################

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_C_INCLUDE_PATH})
set(LIBS ${LIBS} ${MPI_C_LIBRARIES})


################################################################################
# Find PThreads
################################################################################

find_package(Threads REQUIRED)
set(LIBS ${LIBS} ${CMAKE_THREAD_LIBS_INIT})


################################################################################
# Find OpenMP
################################################################################

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


################################################################################
# Find mallocMC
################################################################################

find_package(mallocMC 2.2.0 QUIET)

if(NOT mallocMC_FOUND)
    message(STATUS "Using mallocMC from thirdParty/ directory")
    set(MALLOCMC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../thirdParty/mallocMC")
    find_package(mallocMC 2.2.0 REQUIRED)
endif(NOT mallocMC_FOUND)

include_directories(SYSTEM ${mallocMC_INCLUDE_DIRS})
add_definitions(${mallocMC_DEFINITIONS})
set(LIBS ${LIBS} ${mallocMC_LIBRARIES})


################################################################################
# Find zlib
################################################################################

find_package(ZLIB REQUIRED)
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ZLIB_LIBRARIES})


################################################################################
# Find Boost
################################################################################

find_package(Boost 1.57.0 REQUIRED COMPONENTS program_options regex filesystem
                                              system thread math_tr1
                                              serialization)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})


################################################################################
# PARAM (overwrite) defines
################################################################################

if(PARAM_OVERWRITES)
    foreach(param ${PARAM_OVERWRITES})
       add_definitions(${param})
    endforeach(param)
endif(PARAM_OVERWRITES)


################################################################################
# load cuda_memtest and mpiInfo projects
################################################################################

option(CUDAMEMTEST_ENABLE "Build cuda_memtest and the helper mpiInfo \
                          (allow GPU health test before running PIConGPU)" ON)

if("${PMACC_CUDA_COMPILER}" STREQUAL "clang")
    # cuda_memtest is not compiling with `clang`
    set(CUDAMEMTEST_ENABLE OFF CACHE BOOL "Build cuda_memtest and the helper mpiInfo \
                           (allow GPU health test before running PIConGPU)" FORCE)
    message(STATUS "Disable 'cuda_memtest' build: not supported by clang")
endif()

if(CUDAMEMTEST_ENABLE)
    set(SAME_NVCC_FLAGS_IN_SUBPROJECTS OFF)
    find_path(CUDA_MEMTEST_DIR
            NAMES CMakeLists.txt
            PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../../thirdParty/cuda_memtest"
            DOC "path to cuda_memtest"
            )

    add_subdirectory(${CUDA_MEMTEST_DIR}
                     "${CMAKE_CURRENT_BINARY_DIR}/build_cuda_memtest")


    # mpiInfo utility
    find_path(MPI_INFO_DIR
            NAMES CMakeLists.txt
            PATHS "${CMAKE_CURRENT_SOURCE_DIR}/../mpiInfo"
            DOC "path to mpiInfo"
            )

    add_subdirectory(${MPI_INFO_DIR}
                     "${CMAKE_CURRENT_BINARY_DIR}/build_mpiInfo")
endif()


################################################################################
# PIConGPU options
################################################################################

set(PIC_VERBOSE "1" CACHE STRING
    "Set verbosity level for PIConGPU (default is only physics output)")
add_definitions(-DPIC_VERBOSE_LVL=${PIC_VERBOSE})


################################################################################
# ADIOS
################################################################################

# find adios installation
#   set(ADIOS_USE_STATIC_LIBS ON) # force static linking
find_package(ADIOS 1.10.0)

if(ADIOS_FOUND)
    add_definitions(-DENABLE_ADIOS=1)
    include_directories(SYSTEM ${ADIOS_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${ADIOS_LIBRARIES})
endif(ADIOS_FOUND)


################################################################################
# Warnings
################################################################################

# GNU
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
# ICC
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
# PGI
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minform=inform")
endif()


################################################################################
# libSplash (+ hdf5 due to required headers)
################################################################################

# find libSplash installation
# prefer static libraries over shared ones (but do not force them)
set(Splash_USE_STATIC_LIBS ON)
find_package(Splash 1.6.0 COMPONENTS PARALLEL)

if(Splash_FOUND)
    include_directories(SYSTEM ${Splash_INCLUDE_DIRS})
    list(APPEND Splash_DEFINITIONS "-DENABLE_HDF5=1")
    add_definitions(${Splash_DEFINITIONS})
    set(LIBS ${LIBS} ${Splash_LIBRARIES})
endif(Splash_FOUND)


################################################################################
# PNGwriter
################################################################################

# find PNGwriter installation
find_package(PNGwriter 0.5.6)

if(PNGwriter_FOUND)
    include_directories(SYSTEM ${PNGwriter_INCLUDE_DIRS})
    list(APPEND PNGwriter_DEFINITIONS "-DPIC_ENABLE_PNG=1")
    add_definitions(${PNGwriter_DEFINITIONS})
    set(LIBS ${LIBS} ${PNGwriter_LIBRARIES})
endif(PNGwriter_FOUND)


################################################################################
# ISAAC
################################################################################

find_package(ISAAC 1.3.0 CONFIG QUIET)
if(ISAAC_FOUND)
    message(STATUS "Found ISAAC: ${ISAAC_DIR}")
    SET(ISAAC_STEREO "No" CACHE STRING "Using stereoscopy")
    SET_PROPERTY(CACHE ISAAC_STEREO PROPERTY STRINGS No SideBySide Anaglyph)

    if(${ISAAC_STEREO} STREQUAL "No")
        add_definitions(-DISAAC_STEREO=0)
    endif()
    if(${ISAAC_STEREO} STREQUAL "SideBySide")
        add_definitions(-DISAAC_STEREO=1)
    endif()
    if(${ISAAC_STEREO} STREQUAL "Anaglyph")
        add_definitions(-DISAAC_STEREO=2)
    endif()

    include_directories(SYSTEM ${ISAAC_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${ISAAC_LIBRARIES})

    set(ISAAC_MAX_FUNCTORS "3" CACHE STRING "Max length of the isaac functor chain" )
    set(ISAAC_DEFAULT_WEIGHT "7" CACHE STRING "Default weight of an isaac source" )
    add_definitions(${ISAAC_DEFINITIONS})
    add_definitions(-DISAAC_MAX_FUNCTORS=${ISAAC_MAX_FUNCTORS})
    add_definitions(-DISAAC_FUNCTOR_POW_ENABLED=0)
    add_definitions(-DISAAC_DEFAULT_WEIGHT=${ISAAC_DEFAULT_WEIGHT})

    add_definitions(-DENABLE_ISAAC=1)
else(ISAAC_FOUND)
    if(DEFINED ISAAC_DEPENDENCY_HINTS)
        message(STATUS "ISAAC was found, but detected the following "
                       "problems:" ${ISAAC_DEPENDENCY_HINTS})
    else()
        message(STATUS "Could NOT find ISAAC - set ISAAC_DIR or check your CMAKE_PREFIX_PATH")
    endif()
endif(ISAAC_FOUND)


################################################################################
# Check if PIC_EXTENSION_PATH is relative or absolute
################################################################################

find_path(PIC_EXTENSION_PATH
    NAMES include/simulation_defines/param/components.param
    PATHS "${CMAKE_CURRENT_SOURCE_DIR}"
    DOC "search PATH to extension folder"
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    )

set(PIC_COPY_ON_INSTALL "include/simulation_defines" "submit" CACHE LIST
    "Folder which is copied to installation path during install" )

include_directories(include)
include_directories(BEFORE ${PIC_EXTENSION_PATH}/include)


################################################################################
# Compile & Link PIConGPU
################################################################################

file(GLOB_RECURSE CUDASRCFILES "*.cu")
file(GLOB_RECURSE SRCFILES "*.cpp")

add_library(picongpu-hostonly
    STATIC
    ${SRCFILES}
)

if("${PMACC_CUDA_COMPILER}" STREQUAL "clang")
    add_executable(picongpu
        ${CUDASRCFILES}
    )

    set_target_properties(picongpu PROPERTIES COMPILE_FLAGS ${CLANG_BUILD_FLAGS})
    set_target_properties(picongpu PROPERTIES LINKER_LANGUAGE CXX)
    set_source_files_properties(${CUDASRCFILES} PROPERTIES LANGUAGE CXX)

    target_link_libraries(picongpu ${LIBS} picongpu-hostonly m )
else()
    cuda_add_executable(picongpu
        ${CUDASRCFILES}
        ${SRCFILES}
    )

    target_link_libraries(picongpu ${LIBS} picongpu-hostonly m)
endif()

################################################################################
# Install PIConGPU
################################################################################

install(TARGETS picongpu
         RUNTIME DESTINATION bin)

#file(GLOB scripts_to_copy "${CMAKE_CURRENT_SOURCE_DIR}/scripts/*.sh")
#foreach(f ${scripts_to_copy})
#   GET_FILENAME_COMPONENT(file_name ${f} NAME CACHE)
#   install(FILES "${f}" DESTINATION bin PERMISSIONS OWNER_EXECUTE OWNER_READ
#           OWNER_WRITE GROUP_READ GROUP_EXECUTE)
#endforeach(f)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/scripts/" DESTINATION bin
    FILES_MATCHING PATTERN "*"
    PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_READ GROUP_EXECUTE
    PATTERN .svn EXCLUDE
)

# If the installation prefix does not equal extension path, check if folder must be copied.
# If there is no include folder in installation prefix, also copy all missing folders.
if( (NOT "${CMAKE_INSTALL_PREFIX}" STREQUAL "${PIC_EXTENSION_PATH}") OR
    (NOT EXISTS "${CMAKE_INSTALL_PREFIX}/include"))

    #copy all important subfolders to install folder
    foreach(dir ${PIC_COPY_ON_INSTALL})

      #if source not exists than copy
      if(NOT EXISTS "${CMAKE_INSTALL_PREFIX}/${dir}")
          #copy important subfolders from extension path (default PIConGPU parameter)
          if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/")
            install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/."
              DESTINATION "${CMAKE_INSTALL_PREFIX}/${dir}"
              PATTERN .svn EXCLUDE
              PATTERN .git EXCLUDE
            )
          endif()

          #copy important subfolders from extension path (from extension path)
          if(EXISTS "${PIC_EXTENSION_PATH}/${dir}/")
            install(DIRECTORY "${PIC_EXTENSION_PATH}/${dir}/."
               DESTINATION "${CMAKE_INSTALL_PREFIX}/${dir}"
               PATTERN .svn EXCLUDE
               PATTERN .git EXCLUDE
            )
          endif()
      endif()

    endforeach()
endif()
