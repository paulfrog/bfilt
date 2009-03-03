CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(plane)
ADD_DEFINITIONS(" -O3")

# GSL
SET(BFILT_LIB bfilt)

# Include et Link Directories

IF(APPLE)
  MESSAGE("-- APPLE Configuration")
  INCLUDE_DIRECTORIES(
    /sw/include/
    )
ENDIF(APPLE)

# Executables and "stand-alone " librairies
ADD_EXECUTABLE (plane
  plane.cpp
  example_4.cpp
  )

# Linkage
TARGET_LINK_LIBRARIES(plane
  ${BFILT_LIB})
