# Python bindings for Ultraliser

# Set CMake policy for Python finding (CMP0148)
if(POLICY CMP0148)
    cmake_policy(SET CMP0148 NEW)
endif()

# Enable pybind11's FindPython before finding pybind11
set(PYBIND11_FINDPYTHON ON)

# Find Python first
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Find pybind11
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    find_package(pybind11 QUIET)
endif()

if(NOT pybind11_FOUND)
    message(FATAL_ERROR "pybind11 is required to build the Ultraliser Python bindings")
endif()

set(ULTRALISER_PYTHON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ultraliser_python.cpp
)

pybind11_add_module(ultraliser_python MODULE ${ULTRALISER_PYTHON_SOURCES})

target_link_libraries(ultraliser_python PRIVATE ultraliser)

target_include_directories(ultraliser_python
    PRIVATE
        ${CMAKE_SOURCE_DIR}/ultraliser
        ${CMAKE_BINARY_DIR}/include
)

set_target_properties(ultraliser_python PROPERTIES
    OUTPUT_NAME "_core"
)

install(TARGETS ultraliser_python
        LIBRARY DESTINATION ultraliser
        ARCHIVE DESTINATION ultraliser
        RUNTIME DESTINATION ultraliser)

install(FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ultraliser/__init__.py
        DESTINATION ultraliser)
