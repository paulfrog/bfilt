CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(BFilt)

ADD_DEFINITIONS(" -O3 ")

# GSL
SET(GSL_LIBRARIES gsl)

# ATLAS
SET(LAPACK_LIBS lapack f77blas cblas atlas gfortran)

# Include et Link Directories

IF(APPLE)
  INCLUDE_DIRECTORIES(
    /sw/include/
    ../src
    )
  LINK_DIRECTORIES(
    /sw/lib/  
    )
ELSE(APPLE)
    INCLUDE_DIRECTORIES(
    ../src
    )
ENDIF(APPLE)

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
