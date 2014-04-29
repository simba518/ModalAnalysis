# -*- cmake -*-

# glm: a .obj file loader, which could be download from 
# http://devernay.free.fr/hacks/glm/
# This module defines
# GLM_FOUND, if false, do not try to link to Log4cplus
# GLM_LIBRARIES, include libglm.a, libIL.so and libjpeg.so (the last two are used
# for loading the textures)
# GLM_LIBRARY: libglm.a
# IL_LIBRARY: libIL.so
# JPEG_LIBRARY: libjpeg.so
# GLM_INCLUDE_DIR, where to find log4cplus.hpp

SET(GLM_FOUND "NO")
FIND_PATH(GLM_INCLUDE_DIR glm.h

  /usr/include
  /usr/local/include
  )
IF(NOT GLM_INCLUDE_DIR)
  MESSAGE(FATAL_ERROR "Could not find include file: glm.h")
ENDIF(NOT GLM_INCLUDE_DIR)

FIND_LIBRARY(GLM_LIBRARY

  NAMES libglm.a
  PATHS 
  /usr/lib 
  /usr/local/lib
  )
IF(NOT GLM_LIBRARY)
  MESSAGE(FATAL_ERROR "Could not find GLM library: libglm.a")
ENDIF(NOT GLM_LIBRARY)

FIND_LIBRARY(IL_LIBRARY

  NAMES libIL.so
  PATHS 
  /usr/lib 
  /usr/local/lib
  )
IF(NOT IL_LIBRARY)
  MESSAGE(FATAL_ERROR "Could not find devil library: libIL.so")
ENDIF(NOT IL_LIBRARY)

FIND_LIBRARY(JPEG_LIBRARY

  NAMES libjpeg.so
  PATHS 
  /usr/lib 
  /usr/local/lib
  )
IF(NOT JPEG_LIBRARY)
  MESSAGE(FATAL_ERROR "Could not find jpeg library: libjpeg.so")
ENDIF(NOT JPEG_LIBRARY)

FIND_LIBRARY(PNG_LIBRARY

  NAMES libpng.so
  PATHS 
  /usr/lib
  /usr/local/lib
  /usr/lib/x86_64-linux-gnu
  /usr/local/lib/x86_64-linux-gnu
  )
IF(NOT PNG_LIBRARY)
  MESSAGE(FATAL_ERROR "Could not find jpeg library: libpng.so")
ENDIF(NOT PNG_LIBRARY)

IF (GLM_LIBRARY AND GLM_INCLUDE_DIR AND IL_LIBRARY AND JPEG_LIBRARY AND PNG_LIBRARY)
  SET(GLM_LIBRARIES ${GLM_LIBRARY} ${IL_LIBRARY} ${JPEG_LIBRARY} ${PNG_LIBRARY})
  SET(GLM_FOUND "YES")
ENDIF (GLM_LIBRARY AND GLM_INCLUDE_DIR AND IL_LIBRARY AND JPEG_LIBRARY AND PNG_LIBRARY)

IF (GLM_FOUND)
  MESSAGE(STATUS "Found GLM: ${GLM_LIBRARIES}")
ENDIF (GLM_FOUND)

MARK_AS_ADVANCED(

  GLM_LIBRARIES
  GLM_INCLUDE_DIR
  )
