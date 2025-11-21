####################################################################################################
# Copyright (c) 2016 - 2021
# Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
#
# Author(s)
#       Marwan Abdellah <marwan.abdellah@epfl.ch>
#
# For complete list of authors, please see AUTHORS.md
#
# This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
#
# This library is free software; you can redistribute it and/or modify it under the terms of the
# GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301 USA.
# You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
####################################################################################################

# OpenMP
find_package(OpenMP)

# If OpenMP not found and we're on macOS, try to find Homebrew's libomp
if(NOT OPENMP_FOUND AND APPLE)
    # Try to find Homebrew's libomp
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    
    if(HOMEBREW_LIBOMP_PREFIX)
        message(STATUS "OpenMP not found via standard method, trying Homebrew's libomp at ${HOMEBREW_LIBOMP_PREFIX}")
        
        # Find the OpenMP library
        find_library(OpenMP_C_LIBRARY
            NAMES omp libomp
            PATHS ${HOMEBREW_LIBOMP_PREFIX}/lib
            NO_DEFAULT_PATH
        )
        find_library(OpenMP_CXX_LIBRARY
            NAMES omp libomp
            PATHS ${HOMEBREW_LIBOMP_PREFIX}/lib
            NO_DEFAULT_PATH
        )
        
        if(OpenMP_C_LIBRARY AND OpenMP_CXX_LIBRARY)
            set(OpenMP_CXX_LIBRARIES ${OpenMP_CXX_LIBRARY})
            set(OpenMP_C_LIBRARIES ${OpenMP_C_LIBRARY})
            
            # Set OpenMP compiler flags for Homebrew's libomp
            # -Xpreprocessor -fopenmp is a compiler flag
            # Note: -lomp is NOT included here (it's a linker flag, not a compiler flag)
            if(EXISTS "${HOMEBREW_LIBOMP_PREFIX}/include")
                set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
                set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
            else()
                set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
                set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
            endif()
            
            # Set linker flags separately (these are used during linking, not compilation)
            # The library path is set here, but actual linking is done via link_libraries() below
            set(OpenMP_EXE_LINKER_FLAGS "-L${HOMEBREW_LIBOMP_PREFIX}/lib")
            
            set(OPENMP_FOUND TRUE)
            message(STATUS "Found OpenMP via Homebrew: ${OpenMP_CXX_LIBRARY}")
        endif()
    endif()
endif()

if(OPENMP_FOUND)

    # OpenMP flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

    # Ultaliser OpenMP
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DULTRALISER_USE_OPENMP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_USE_OPENMP")
    set(ULTRALISER_USE_OPENMP TRUE)
    message(STATUS "Found OpenMP " ${OpenMP_C_FLAGS} ", "
                                   ${OpenMP_CXX_FLAGS}, ", "
                                   ${OpenMP_CXX_LIBRARIES})
    # Link libraries
    link_libraries(${OpenMP_CXX_LIBRARIES})
else(OPENMP_FOUND)
    message(STATUS "OpenMP NOT Found")
    set(ULTRALISER_USE_OPENMP FALSE)
endif(OPENMP_FOUND)
