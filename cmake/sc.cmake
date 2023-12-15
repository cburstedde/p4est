# provides imported target SC::SC
include(FetchContent)

# checking if we are building from a git repository or using an archive file
# if using an archive, we should be able to build p4est without git installed.
if(EXISTS "${PROJECT_SOURCE_DIR}/.git")

  # make sure libsc git submodule is up to date
  include(${CMAKE_CURRENT_LIST_DIR}/GitSubmodule.cmake)
  git_submodule("${PROJECT_SOURCE_DIR}/sc")

endif()

# check that libsc source directory is populated
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/sc/CMakeLists.txt")
  message(FATAL_ERROR "Something is wrong, libsc sources directory is not valid. Please check your archive file.")
endif()

FetchContent_Declare(SC SOURCE_DIR ${PROJECT_SOURCE_DIR}/sc)

FetchContent_MakeAvailable(SC)
