CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(ar_process)
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
ADD_EXECUTABLE (ar_process
  ar_process.cpp
  example_1.cpp
  )

# Linkage
TARGET_LINK_LIBRARIES(ar_process
  ${BFILT_LIB})
