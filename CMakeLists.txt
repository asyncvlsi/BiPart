cmake_minimum_required(VERSION 3.16)

project(bipart
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

message(STATUS "Detecting environment variable $ENV{ACT_HOME}...")
if (DEFINED ENV{ACT_HOME})
    message(STATUS "Environment variable ACT_HOME detected: " $ENV{ACT_HOME})
else ()
    message(FATAL_ERROR "Environment variable $ENV{ACT_HOME} not found")
endif ()
include_directories($ENV{ACT_HOME}/include)
link_directories($ENV{ACT_HOME}/lib)
find_library(PHYDB_LIBRARY
    NAMES phydb
    PATHS $ENV{ACT_HOME}/lib
    REQUIRED)
if (NOT PHYDB_LIBRARY)
    link_directories($ENV{ACT_HOME}/lib)
    set(PHYDB_LIBRARY $ENV{ACT_HOME}/lib)
    #message(FATAL_ERROR "libphydb.a is not found")
else ()
  message(STATUS "libphydb.a found")
endif()

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

message(STATUS "Detecting environment variable LEF_ROOT...")
set(LEF_ROOT "")
if (DEFINED ENV{LEF_ROOT})
    message(STATUS "Environment variable LEF_ROOT detected: " $ENV{LEF_ROOT})
    include_directories($ENV{LEF_ROOT}/include)
    link_directories($ENV{LEF_ROOT}/lib)
    set(LEF_ROOT $ENV{LEF_ROOT})
else ()
    message(STATUS "Environment variable LEF_ROOT not found, using ACT_HOME instead")
    include_directories($ENV{ACT_HOME}/include)
    link_directories($ENV{ACT_HOME}/lib)
    set(LEF_ROOT $ENV{ACT_HOME})
endif ()
#include_directories($ENV{LEF_ROOT}/include)
#link_directories($ENV{LEF_ROOT}/lib)
find_library(LEF_LIBRARY
    NAMES lef
    PATHS $ENV{LEF_ROOT}/lib
    REQUIRED)
if (NOT LEF_LIBRARY)
    message(FATAL_ERROR "liblef.a is not found")
else ()
    message(STATUS "liblef.a found")
endif()

message(STATUS "Detecting environment variable DEF_ROOT...")
set(DEF_ROOT "")
if (DEFINED ENV{DEF_ROOT})
    message(STATUS "Environment variable DEF_ROOT detected: " $ENV{DEF_ROOT})
    include_directories($ENV{DEF_ROOT}/include)
    link_directories($ENV{DEF_ROOT}/lib)
    set(DEF_ROOT $ENV{DEF_ROOT})
else ()
    message(STATUS "Environment variable DEF_ROOT not found, using ACT_HOME instead")
    include_directories($ENV{ACT_HOME}/include)
    link_directories($ENV{ACT_HOME}/lib)
    set(DEF_ROOT $ENV{ACT_HOME})
endif ()

find_library(DEF_LIBRARY
    NAMES def
    PATHS $ENV{DEF_ROOT}/lib
    REQUIRED)
if (NOT DEF_LIBRARY)
    message(FATAL_ERROR "libdef.a is not found")
else ()
    message(STATUS "libdef.a found")
endif()



message(STATUS "Detecting environment variable $ENV{GALOIS_INCLUDE}...")
if (DEFINED ENV{GALOIS_INCLUDE})
    message(STATUS "Environment variable GALOIS_INCLUDE detected: " $ENV{GALOIS_INCLUDE})
else ()
    message(FATAL_ERROR "Environment variable $ENV{GALOIS_INCLUDE} not found")
endif ()



message(STATUS "Detecting Galois libraries...")
set(Galois_FOUND TRUE)
find_library(
    GALOIS_EDA_LIBRARY
    NAMES galois_eda
    PATHS $ENV{ACT_HOME}/lib
)
if (GALOIS_EDA_LIBRARY)
    set(Galois_LIBRARIES ${GALOIS_EDA_LIBRARY})
    message(STATUS "Found libgalois_eda.a: " ${GALOIS_EDA_LIBRARY})
else()
    message(STATUS "Cannot find libgalois_eda.a")
    set(Galois_FOUND FALSE)
endif()
find_library(
    GALOIS_SHMEM_LIBRARY
    NAMES galois_shmem
    PATHS $ENV{ACT_HOME}/lib
)
if (GALOIS_SHMEM_LIBRARY)
    set(Galois_LIBRARIES ${Galois_LIBRARIES} ${GALOIS_SHMEM_LIBRARY})
    message(STATUS "Found libgalois_shmem.a: " ${GALOIS_SHMEM_LIBRARY})
else()
    message(STATUS "Cannot find libgalois_shmem.a")
    set(Galois_FOUND FALSE)
endif()
if (UNIX AND (NOT APPLE))
    find_library(
        NUMA_LIBRARY
        NAMES numa
    )
    if (NUMA_LIBRARY)
        set(Galois_LIBRARIES ${Galois_LIBRARIES} ${NUMA_LIBRARY})
        message(STATUS "Found libnuma: " ${NUMA_LIBRARY})
    else()
        message(STATUS "Cannot find libnuma")
        set(Galois_FOUND FALSE)
    endif()
endif()



include_directories($ENV{GALOIS_INCLUDE})
#include_directories($ENV{PHYDB_INCLUDE})
#include_directories("/net/ohm/export/iss/sepideh/act/install/include")
#include_directories("/net/ohm/export/cdgc/michael/GaloisCpp-github/Galois/lonestar/liblonestar/include")
#include_directories("/net/ohm/export/iss/sepideh/Galois/libgalois/include")
#include_directories("/net/ohm/export/iss/sepideh/Galois/build/libgalois/include")
#include_directories("/net/ohm/export/iss/sepideh/Galois/lonestar/liblonestar/include")
#include_directories("/org/centers/cdgc/sw/llvm-10.0.0/include")
link_directories($ENV{GALOIS_LIB})
find_library(GALOIS_LIB
    NAMES galois_shmem
    PATHS $ENV{GALOIS_LIB}
    REQUIRED)
if (NOT GALOIS_LIB)
    message(FATAL_ERROR "libgalois_shmem.a is not found")
else ()
    message(STATUS "libgalois_shmem.a found")
endif()



find_package(Boost 1.69.0 REQUIRED)
message(STATUS "Boost library path: ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost include path: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost libs: ${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS})

# Set a default build type if none was specified
set(default_build_type "RELEASE")
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}")
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "DEBUG" "RELEASE")
endif ()

#add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic)
#set(CMAKE_CXX_FLAGS_DEBUG "-g -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-function -Wno-shadow -Wno-unused-parameter")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-function -Wno-shadow -Wno-unused-parameter")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG -Wno-unused-variable -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-function -Wno-shadow -Wno-unused-parameter")

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Set the output directory of static libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# Set the output directory of executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

include_directories(include)
FILE(GLOB SOURCES "src/*.cpp")
#list(FILTER CPP_SOURCES EXCLUDE REGEX main)
#FILE(GLOB C_SOURCES "src/*.c")
add_subdirectory(test)

add_library(bipart 
    STATIC
    ${SOURCES})


target_link_libraries(bipart
  lef def phydb ${Galois_LIBRARIES} pthread)

#add_executable(bipart_test
#  test/main.cpp)

#target_link_libraries(bipart_test
#  lef def galois_shmem galois_eda numa bipart phydb pthread)




######################################
# specify the installation directory: ${ACT_HOME}
######################################

message(STATUS "Detecting environment variable ACT_HOME...")
if (DEFINED ENV{ACT_HOME})
    message(STATUS "Environment variable ACT_HOME detected: " $ENV{ACT_HOME})
else ()
    message(FATAL_ERROR "Environment variable ACT_HOME not found")
endif ()
message(STATUS "Changing Installation directory to ACT_HOME")
set(CMAKE_INSTALL_PREFIX $ENV{ACT_HOME} CACHE PATH "installation path" FORCE)
message(STATUS "Current installation directory: " ${CMAKE_INSTALL_PREFIX})


######################################
# Install header files
######################################
install(DIRECTORY include/ DESTINATION include/bipart COMPONENT Development FILES_MATCHING PATTERN "*.h")

######################################
# Install library
######################################
install(TARGETS bipart
    DESTINATION lib)
