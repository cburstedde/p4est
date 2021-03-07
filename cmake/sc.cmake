# provides imported target SC::SC
include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)

set(sc_external true CACHE BOOL "build sc library" FORCE)

git_submodule("${PROJECT_SOURCE_DIR}/sc")
