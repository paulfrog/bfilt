CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(Ornstein)
ADD_DEFINITIONS(" -O3")

# GSL
SET(BFILT_LIB bfilt)

# Include et Link Directories

IF(APPLE)
  MESSAGE("-- Apple Configuration")
  INCLUDE_DIRECTORIES(
    /sw/include/
    )
ENDIF(APPLE)

# Executables and "stand-alone " librairies
ADD_EXECUTABLE (Ornstein
  ornstein_uhlenbeck.cpp
  example_2.cpp
  )

# Linkage
TARGET_LINK_LIBRARIES(Ornstein
  ${BFILT_LIB}     
  )
