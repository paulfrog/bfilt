CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(BFilt)
ADD_DEFINITIONS(" -O3 ")

# GSL
SET(GSL_LIBRARIES gsl)

# ATLAS
SET(LAPACK_LIBS lapack f77blas cblas atlas gfortran)

# Include et Link Directories

INCLUDE_DIRECTORIES(
  /sw/include/
  /usr/local/include
  ../src
  )
LINK_DIRECTORIES(
  /sw/lib/  
  /usr/local/lib 
  )


# Executables and "stand-alone " librairies
ADD_LIBRARY (bfilt SHARED
  simulator.cpp
  filter_tools.cpp
  filter.cpp
  extended_kalman_filter.cpp
  gaussian_model.cpp
  unscented_kalman_filter.cpp
  sisr_filter.cpp
  local_linearization_filter.cpp
  thgl_filter.cpp
  )


# Linkage
TARGET_LINK_LIBRARIES(bfilt
  ${GSL_LIBRARIES}     
  ${LAPACK_LIBS})


# INSTALLATION
INSTALL(TARGETS bfilt DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL(DIRECTORY ../src/bfilt DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
  FILES_MATCHING PATTERN "*.h")

#documentation

##Generation of the documentation using doxygen
FIND_PATH(DOXYGEN_DIR doxygen
  /usr/local/bin
  /usr/bin
  )
IF (DOXYGEN_DIR)
  MESSAGE(STATUS "Doxygen generation of the documentation")
  EXEC_PROGRAM(
    "cd ../doc ; ${DOXYGEN_DIR}/doxygen Doxyfile"
    OUTPUT_VARIABLE MY_OUTPUT
    )
ELSE (DOXYGEN_DIR)
  MESSAGE(STATUS "Doxygen not found, I cannot generate the documentation")
ENDIF (DOXYGEN_DIR)
