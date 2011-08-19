PROJECT(minirosetta)

# if you don't want the full compiler output, remove the following line
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# this seems to fix some linking errors
SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cq <TARGET> <LINK_FLAGS> <OBJECTS>")
SET(CMAKE_CXX_ARCHIVE_APPEND "<CMAKE_AR> q  <TARGET> <LINK_FLAGS> <OBJECTS>")

# include / link dirs
INCLUDE_DIRECTORIES(../..)
INCLUDE_DIRECTORIES(../../src)
INCLUDE_DIRECTORIES(../../external/cxxtest)

# external libraries
CMAKE_POLICY( SET CMP0015 NEW )
INCLUDE_DIRECTORIES(../../external/boost_1_46_1)
LINK_DIRECTORIES(../../external/boost_1_46_1)

INCLUDE_DIRECTORIES(../../external/include)
LINK_DIRECTORIES(../../external/lib)

# Platform-specific includes
if(APPLE)
	INCLUDE_DIRECTORIES(../../src/platform/macos)
	ADD_DEFINITIONS(-DMAC)
ELSEIF(UNIX)
	INCLUDE_DIRECTORIES(../../src/platform/linux)
	ADD_DEFINITIONS(-DLINUX)
ELSEIF(WIN32)
	INCLUDE_DIRECTORIES(../../src/platform/windows)
	ADD_DEFINITIONS(-DWIN32)
ENDIF(APPLE)
