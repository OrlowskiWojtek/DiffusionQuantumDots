cmake_minimum_required(VERSION 3.10)

# Create the library target
add_library(WalkersVisualiser
    src/visualiser.cpp
)

# Use C++20
target_compile_features(WalkersVisualiser PUBLIC cxx_std_20)

# Platform-specific compile options
if (APPLE)
    target_compile_options(WalkersVisualiser PRIVATE -Wall -Wfatal-errors -g -O3 -DGL_SILENCE_DEPRECATION)
elseif (WIN32)
    target_compile_options(WalkersVisualiser PRIVATE -DNOMINMAX /EHsc)
else()
    target_compile_options(WalkersVisualiser PRIVATE -Wall -Wfatal-errors -g -O3 -Wno-unused-result -Wno-unknown-pragmas)
endif()

# OpenMP support
find_package(OpenMP 3 QUIET)
if(OpenMP_FOUND)
    target_link_libraries(WalkersVisualiser PUBLIC OpenMP::OpenMP_CXX)
endif()

# morphologica fonts
target_compile_definitions(WalkersVisualiser PUBLIC
    MORPH_FONTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/morphologica/fonts"
    FLT=float
    COMPILE_PLOTTING
)

# morphologica path setup
set(MORPH_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/morphologica" CACHE PATH "The path to morphologica")

# Required libraries
set(OpenGL_GL_PREFERENCE GLVND)

find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED)
find_package(Freetype REQUIRED)
find_package(rapidxml QUIET)

# Fallback for rapidxml
if(NOT rapidxml_FOUND)
    message("-- Using the bundled version of rapidxml")
    target_include_directories(WalkersVisualiser PRIVATE "${MORPH_INCLUDE_PATH}/include/rapidxml-1.13")
endif()

# Include paths
target_include_directories(WalkersVisualiser
    PUBLIC
        ${MORPH_INCLUDE_PATH}/include
        ${MORPH_INCLUDE_PATH}
        ${ARMADILLO_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${OPENGL_INCLUDE_DIR}
        ${GLFW3_INCLUDE_DIR}
        ${FREETYPE_INCLUDE_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link libraries
target_link_libraries(WalkersVisualiser
    PUBLIC
        ${ARMADILLO_LIBRARIES}
        ${HDF5_C_LIBRARIES}
        OpenGL::GL
        Freetype::Freetype
        glfw
)
