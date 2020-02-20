CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(cam7_kessler)
ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{NETCDF})

#-----------------------------------------------------------------------------
#
# Create list of SCHEME_FILES, HOST_FILES, and SUITE_FILES
# Paths should be relative to CMAKE_SOURCE_DIR (this file's directory)
#
#-----------------------------------------------------------------------------

#NOTE -- FILES WITH DDTs MUST APPEAR BEFORE THEY ARE USED

LIST(APPEND DDT_FILES      "../src/cam7_ddt_files.txt" )
LIST (APPEND SCHEME_FILES  ${DDT_FILES} "simple_schemes.txt" )
LIST(APPEND HOST_FILES     "../src/constituents" "../src/physconst" "../src/ppgrid" )
LIST(APPEND SUITE_FILES    "../../atmos_phys/suite_kessler.xml")
LIST(APPEND UTILITY_FILES  "../../atmos_phys/kessler/geopotential_t.F90" )

# HOST is the name of the executable we will build.
SET(HOST "${CMAKE_PROJECT_NAME}")
# HOST_PROGRAM is the name of the program source files (Fortran and metadata)
SET(HOST_PROGRAM "../src/cam7_kessler")

#-----------------------------------------------------------------------------
#
# Set where the CCPP Framework lives
#
#-----------------------------------------------------------------------------
get_filename_component(kessler_ROOT "${CMAKE_SOURCE_DIR}" DIRECTORY)
set(CCPP_ROOT "${kessler_ROOT}/../ccpp-framework")

#-----------------------------------------------------------------------------
##############################################################################
#
# End of project-specific input
#
##############################################################################
#-----------------------------------------------------------------------------

# By default, no verbose output
SET(VERBOSITY 0 CACHE STRING "Verbosity level of output (default: 0)")
# By default, generated caps go in ccpp subdir
SET(CCPP_CAP_FILES "${CMAKE_BINARY_DIR}/ccpp" CACHE
  STRING "Location of CCPP-generated cap files")

SET(CCPP_FRAMEWORK ${CCPP_ROOT}/scripts)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    #message(STATUS "Setting build type to 'Debug' as none was specified.")
    #set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                 "MinSizeRel" "RelWithDebInfo")
endif()

ADD_COMPILE_OPTIONS(-O0)

ADD_EXECUTABLE(${HOST} ${HOST_PROGRAM}.F90)
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
#  gfortran
#  MESSAGE("gfortran being used.")
  ADD_COMPILE_OPTIONS(-fcheck=all)
  ADD_COMPILE_OPTIONS(-fbacktrace)
  ADD_COMPILE_OPTIONS(-ffpe-trap=zero)
  ADD_COMPILE_OPTIONS(-finit-real=nan)
  ADD_COMPILE_OPTIONS(-ggdb)
  ADD_COMPILE_OPTIONS(-ffree-line-length-none)
  ADD_COMPILE_OPTIONS(-cpp)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
#  ifort
#  MESSAGE("ifort being used.")
  #ADD_COMPILE_OPTIONS(-check all)
  ADD_COMPILE_OPTIONS(-fpe0)
  ADD_COMPILE_OPTIONS(-warn)
  ADD_COMPILE_OPTIONS(-traceback)
  ADD_COMPILE_OPTIONS(-debug extended)
  ADD_COMPILE_OPTIONS(-fpp)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
#  pgf90
#  MESSAGE("pgf90 being used.")
  ADD_COMPILE_OPTIONS(-g)
  ADD_COMPILE_OPTIONS(-Mipa=noconst)
  ADD_COMPILE_OPTIONS(-traceback)
  ADD_COMPILE_OPTIONS(-Mfree)
  ADD_COMPILE_OPTIONS(-Mfptrap)
  ADD_COMPILE_OPTIONS(-Mpreprocess)
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#-----------------------------------------------------------------------------
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CCPP_FRAMEWORK}/cmake")
#-----------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
else(OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif()

# Create metadata and source file lists
#FOREACH (FILE ${UTILITY_FILES})
#  FILE(STRINGS ${FILE} U_FILENAMES)
#  LIST(APPEND UTILITY_FILENAMES ${U_FILENAMES})
#ENDFOREACH(FILE)

FOREACH (FILE ${SCHEME_FILES})
  FILE(STRINGS ${FILE} FILENAMES)
  LIST(APPEND SCHEME_FILENAMES ${FILENAMES})
ENDFOREACH(FILE)
string(REPLACE ";" "," SCHEME_METADATA "${SCHEME_FILES}")

FOREACH(FILE ${SCHEME_FILENAMES})
  # target_sources prefers absolute pathnames
  string(REPLACE ".meta" ".F90" TEMP "${FILE}")
  get_filename_component(ABS_PATH "${TEMP}" ABSOLUTE)
  MESSAGE(STATUS "HOST " ${HOST} " ABS_PATH " ${ABS_PATH})
  target_sources(${HOST} PUBLIC ${ABS_PATH})
ENDFOREACH(FILE)

FOREACH (FILE ${DDT_FILES})
  FILE(STRINGS ${FILE} FILENAMES)
  LIST(APPEND DDT_FILENAMES ${FILENAMES})
ENDFOREACH(FILE)
string(REPLACE ".meta" ";" DDT_FILENAMES "${DDT_FILENAMES}")

LIST(APPEND DDT_HOST_FILES ${DDT_FILENAMES} ${HOST_FILES})

FOREACH(FILE ${DDT_HOST_FILES})
  LIST(APPEND DDT_HOST_METADATA "${FILE}.meta")
  # target_sources prefers absolute pathnames
  get_filename_component(ABS_PATH "${FILE}.F90" ABSOLUTE)
  LIST(APPEND HOST_SOURCE   "${ABS_PATH}")
ENDFOREACH(FILE)
target_sources(${HOST} PUBLIC ${HOST_SOURCE})
string(REPLACE ";" ".meta," DDT_HOST_METADATA "${DDT_HOST_FILES}")
set(DDT_HOST_METADATA "${DDT_HOST_METADATA}.meta,${HOST_PROGRAM}.meta")

string(REPLACE ";" "," SUITE_XML "${SUITE_FILES}")

# Run ccpp_capgen
set(CAPGEN_CMD "${CCPP_FRAMEWORK}/ccpp_capgen.py")
list(APPEND CAPGEN_CMD "--host-files")
list(APPEND CAPGEN_CMD "${DDT_HOST_METADATA}")
list(APPEND CAPGEN_CMD "--scheme-files")
list(APPEND CAPGEN_CMD "${SCHEME_METADATA}")
list(APPEND CAPGEN_CMD "--suites")
list(APPEND CAPGEN_CMD "${SUITE_XML}")
list(APPEND CAPGEN_CMD "--output-root")
list(APPEND CAPGEN_CMD "${CCPP_CAP_FILES}")
list(APPEND CAPGEN_CMD "--generate-host-cap")
list(APPEND CAPGEN_CMD "--host-name")
list(APPEND CAPGEN_CMD "${HOST}")
  list(APPEND CAPGEN_CMD "--verbose")
while (VERBOSITY GREATER 0)
  list(APPEND CAPGEN_CMD "--verbose")
  MATH(EXPR VERBOSITY "${VERBOSITY} - 1")
endwhile ()
string(REPLACE ";" " " CAPGEN_STRING "${CAPGEN_CMD}")
MESSAGE(STATUS "Running: ${CAPGEN_STRING}")
EXECUTE_PROCESS(COMMAND ${CAPGEN_CMD} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE CAPGEN_OUT ERROR_VARIABLE CAPGEN_OUT RESULT_VARIABLE RES)
MESSAGE(STATUS "${CAPGEN_OUT}")
if (RES EQUAL 0)
  MESSAGE(STATUS "CCPP cap generation completed")
else(RES EQUAL 0)
  MESSAGE(FATAL_ERROR "CCPP cap generation FAILED: result = ${RES}")
endif(RES EQUAL 0)

# Read the list of files from capfiles.txt and set to CCPP_CAPS
file(STRINGS ${CCPP_CAP_FILES}/capfiles.txt CCPP_CAPS)
target_sources(${HOST} PUBLIC ${CCPP_CAPS} ${UTILITY_FILENAMES})

SET(NETCDF_F90 "YES")
FIND_PACKAGE(NetCDF)

INCLUDE_DIRECTORIES(${NETCDF_INCLUDES})
INCLUDE_DIRECTORIES(${CCPP_CAP_FILES})

TARGET_LINK_LIBRARIES(cam7_kessler ${NETCDF_LIBRARIES})

set_target_properties(${HOST} PROPERTIES
                              COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                              LINK_FLAGS "${CMAKE_Fortran_FLAGS}")
