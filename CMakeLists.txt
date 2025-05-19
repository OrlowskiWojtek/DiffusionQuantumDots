cmake_minimum_required(VERSION 3.22)

project(DiffusionQuantumMC VERSION 1.0)

# specify the C++ standard
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_23)

# List of compilers
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")

# Adding compiler flags
target_compile_options(compiler_flags INTERFACE
                        "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wformat=2;-Wunused;-Werror>>"
                        "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>")

configure_file(DQMCConfig.h.in DQMCConfig.h)

add_subdirectory(src/UnitHandler)
add_subdirectory(src/DiffusionSolver)
add_subdirectory(src/FileLoader)

## BOOST setup
FIND_PACKAGE( Boost 1.74 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

if(BOOST_FOUND)
    message("Found boost, using boost ver " ${Boost_VERSION_STRING})
endif()

## Parallel TBB setup
find_package(TBB REQUIRED)

add_compile_options(-ltbb -fopenmp)
add_definitions(-DPSTL_USE_PARALLEL_POLICIES)

## Armadillo setup
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

## openMP setup

add_executable(${PROJECT_NAME} src/main.cpp)

target_link_libraries(DiffusionQuantumMC DiffusionSolver FileLoader compiler_flags TBB::tbb ${ARMADILLO_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR})
