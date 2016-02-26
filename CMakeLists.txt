#
# This is a CMake makefile.  You can find the cmake utility and
# information about it at http://www.cmake.org
#


#-------------------------------------------------------------
# set default values of some CMAKE variables on our project
#-------------------------------------------------------------


#EXEC_PROGRAM(whoami OUTPUT_VARIABLE BINPATH)
SET(BIN_PATH "./bin/")

SET(3DIMSEG ./header/ )

set (CMAKE_BUILD_TYPE Release
	CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." )


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set (CMAKE_CXX_FLAGS "-Wno-deprecated -Wno-bool-conversions -fopenmp" 
		CACHE STRING "Flags used by the compiler during all build types." )	

else ()
	set (CMAKE_CXX_FLAGS "-fopenmp -Wno-deprecated" 
		CACHE STRING "Flags used by the compiler during all build types." )	


endif()

set (CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -pedantic -DFPE_ABORT"
	CACHE STRING "Flags used by the compiler during C++ debug builds." )


set (EXECUTABLE_OUTPUT_PATH ${BIN_PATH} 
	CACHE PATH "Single output directory for building all executables.")






#-------------------------------------------------------------
#					starting project
#-------------------------------------------------------------

PROJECT (3dimagetoolkit CXX)




INCLUDE(CMakeDependentOption)

cmake_minimum_required(VERSION 2.8)




#------------------------------------
# looking for STL Library
#------------------------------------

find_library (LIB_DL NAMES dl)
find_library (LIB_M NAMES m)

link_directories (${LIB_DL//libdl.so/})
link_directories (${LIB_M//libm.so/})




#------------------------------------
# looking for VTK libraries
#------------------------------------

FIND_PACKAGE(VTK REQUIRED)
MESSAGE("Using VTK version: ${VTK_VERSION}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVTK_VERSION_MAJOR=${VTK_VERSION_MAJOR}") 
IF(VTK_VERSION_MAJOR VERSION_LESS 6)
    IF(NOT VTK_USE_RENDERING)
        MESSAGE(FATAL_ERROR "Example ${PROJECT_NAME} requires VTK_USE_RENDERING.")
    ENDIF(NOT VTK_USE_RENDERING)
ENDIF()
INCLUDE(${VTK_USE_FILE})



#------------------------------------
# looking for fftw libraries
#------------------------------------
# Find the native FFTW includes and library
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.
#  FFTW_PARALLEL

if (FFTW_INCLUDES)
	# Already in cache, be silent
	set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h)

find_library (FFTW_LIBRARIES NAMES fftw3)

# find fftw_omp_library
find_library (FFTW_PARALLEL_LIB NAMES fftw3_omp)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_PARALLEL_LIB FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES FFTW_PARALLEL_LIB LIB_DL LIB_M VTK_DIR)






INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples/data DESTINATION ${EXECUTABLE_OUTPUT_PATH})


#------------------------------------
# create EXECUTABLES
#------------------------------------

OPTION(BUILD_EXECUTABLES 
		"Build all example executables, set OFF if you want to build only some executable" 
		ON)

CMAKE_DEPENDENT_OPTION ( BUILD_SHOW_IMAGE 
		"Build an executable useful to show a 3d image containing a level set function" 
		OFF
		"NOT BUILD_EXECUTABLES" 
		ON	)

CMAKE_DEPENDENT_OPTION ( BUILD_RSFE_SPLITBREGMAN 
		"Build an executable useful to test rsfe algorithm on a 3d image" 
		OFF
		"NOT BUILD_EXECUTABLES" 
		ON	)

CMAKE_DEPENDENT_OPTION ( BUILD_IMAGE_TOOLKIT
		"Build an executable useful to use all members of the library" 
		OFF
		"NOT BUILD_EXECUTABLES" 
		ON	)

CMAKE_DEPENDENT_OPTION ( BUILD_MYEXEC
		"Build an executable customizable from user" 
		OFF
		"NOT BUILD_EXECUTABLES" 
		ON	)

CMAKE_DEPENDENT_OPTION ( BUILD_TEST
		"Build an executable to test the library" 
		OFF
		"NOT BUILD_EXECUTABLES" 
		OFF	)

# 1
if(BUILD_SHOW_IMAGE)
	ADD_EXECUTABLE(show_image ./examples/showing.cxx)

	# link EXECUTABLE to our library
	TARGET_LINK_LIBRARIES(show_image ${VTK_LIBRARIES} )


endif()


# 2
if(BUILD_RSFE_SPLITBREGMAN)
	ADD_EXECUTABLE(rsfe_splitbregman ./examples/rsfe_splitbregman.cxx)

	# link EXECUTABLE to our library
	TARGET_LINK_LIBRARIES(rsfe_splitbregman
							${VTK_LIBRARIES} ${FFTW_PARALLEL_LIB} ${FFTW_LIBRARIES} )


endif()


# 3
if(BUILD_IMAGE_TOOLKIT)
	ADD_EXECUTABLE(image_toolkit ./examples/image_toolkit.cxx)

	# link EXECUTABLE to our library
	TARGET_LINK_LIBRARIES(image_toolkit
							${VTK_LIBRARIES} ${FFTW_PARALLEL_LIB} ${FFTW_LIBRARIES} )

endif()

# 4
if(BUILD_MYEXEC)
	ADD_EXECUTABLE(myexec ./examples/myexec.cxx)

	# link EXECUTABLE to our library
	TARGET_LINK_LIBRARIES(myexec
							${VTK_LIBRARIES} ${FFTW_PARALLEL_LIB} ${FFTW_LIBRARIES} )

endif()

# 5
if(BUILD_TEST)
	ADD_EXECUTABLE(test ./examples/test.cxx)

	# link EXECUTABLE to our library
	TARGET_LINK_LIBRARIES(test
							${VTK_LIBRARIES} ${FFTW_PARALLEL_LIB} ${FFTW_LIBRARIES} )

endif()





#------------------------------------------------------------------
#						make doc
#------------------------------------------------------------------

#-- Add an Option to toggle the generation of the API documentation
option(BUILD_DOCUMENTATION 
	"Use Doxygen to create the HTML based API documentation typing make Docs on the shell" 
	OFF)

if(BUILD_DOCUMENTATION)
	FIND_PACKAGE(Doxygen)
	if (NOT DOXYGEN_FOUND)
	message(FATAL_ERROR 
		"Doxygen is needed to build the documentation. Please install it correctly")
	endif()
	#-- Configure the Template Doxyfile for our specific project
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
	#-- Add a custom target to run Doxygen when ever the project is built
	add_custom_target (doc
	                    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
	                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	                    COMMENT "Generating API documentation with Doxygen" VERBATIM)
	# IF you do NOT want the documentation to be generated EVERY time you build the project
	# then leave out the 'ALL' keyword from the above command.
endif()



#--------------------------------------------------------------------
#        				make cmakereset
#--------------------------------------------------------------------

add_custom_target (cmakereset
			COMMAND rm -R CMakeCache.txt CMakeFiles/ cmake_install.cmake
			)


#---------------------------------------------------------------------
#						make install
#---------------------------------------------------------------------


FILE(GLOB header 
		"${CMAKE_CURRENT_SOURCE_DIR}/header/*.hxx")

INSTALL(FILES ${header} DESTINATION include)



#INSTALL(TARGETS image_toolkit rsfe_splitbregman show_image
#  RUNTIME DESTINATION bin)




