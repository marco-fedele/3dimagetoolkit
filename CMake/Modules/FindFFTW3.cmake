# - Try to find FFTW3
# Once done this will define
#  FFTW3_FOUND - System has fftw3
#  FFTW3_INCLUDE_DIRS - The fftw3 include directories
#  FFTW3_LIBRARIES - The libraries needed to use fftw3
#  FFTW3_DEFINITIONS - Compiler switches required for using fftw3

find_package(PkgConfig)
pkg_check_modules(PC_LIBXML QUIET fftw3)
set(FFTW3_DEFINITIONS ${PC_FFTW3_CFLAGS_OTHER})

set( FFTW3_ROOT $ENV{FFTW3_ROOT} )
#message("found FFTW3_ROOT from env $ENV{FFTW3_ROOT} and set to ${FFTW3_ROOT}" )

if(FFTW3_ROOT)
    find_path(FFTW3_INCLUDE_DIR fftw3.h 
          HINTS ${FFTW3_ROOT}
		/opt/local/
		/usr/local/
          PATH_SUFFIXES include)
else()
    message("Please specify FFTW3_ROOT variable as the root of the fftw3 installation")
endif()
#message("found FFTW3 includes at ${FFTW3_INCLUDE_DIR}" )

find_library(FFTW3_LIBRARY fftw3
             HINTS ${FFTW3_ROOT}
		/opt/local/
		/usr/local/
             PATH_SUFFIXES lib
            )

find_library(FFTW3_PARALLEL_LIB fftw3_omp 
             HINTS ${FFTW3_ROOT}
		/opt/local/
		/usr/local/
             PATH_SUFFIXES lib
            )

set(FFTW3_LIBRARIES ${FFTW3_LIBRARY} ${FFTW3_PARALLEL_LIB})
set(FFTW3_INCLUDE_DIR ${FFTW3_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(fftw3  DEFAULT_MSG
                                  FFTW3_LIBRARY FFTW3_PARALLEL_LIB FFTW3_INCLUDE_DIR)

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY FFTW_PARALLEL_LIB)
