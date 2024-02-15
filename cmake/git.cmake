# --- extract version from Git

set(PROJECT_MAJOR 0)
set(PROJECT_MINOR 0)
set(PROJECT_PATCH 0)
set(PROJECT_VERSION 0.0.0)
find_program(GIT_VERSION_GEN NAMES git-version-gen
             PATHS ${PROJECT_SOURCE_DIR}/build-aux NO_DEFAULT_PATH)
if(GIT_VERSION_GEN)
  execute_process(COMMAND ${GIT_VERSION_GEN} .tarball-version
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE _err
    OUTPUT_VARIABLE git_version
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
if(_err EQUAL 0)
  if(git_version MATCHES
                 "^(0|[1-9][0-9]*)[.](0|[1-9][0-9]*)[.](0|[1-9][0-9]*)[.].*")
    set(PROJECT_MAJOR "${CMAKE_MATCH_1}")
    set(PROJECT_MINOR "${CMAKE_MATCH_2}")
    set(PROJECT_PATCH "${CMAKE_MATCH_3}")
    set(PROJECT_VERSION ${PROJECT_MAJOR}.${PROJECT_MINOR}.${PROJECT_PATCH}.999)
  elseif(git_version MATCHES
                 "^(0|[1-9][0-9]*)[.](0|[1-9][0-9]*)[.](0|[1-9][0-9]*)")
    set(PROJECT_MAJOR "${CMAKE_MATCH_1}")
    set(PROJECT_MINOR "${CMAKE_MATCH_2}")
    set(PROJECT_PATCH "${CMAKE_MATCH_3}")
    set(PROJECT_VERSION ${PROJECT_MAJOR}.${PROJECT_MINOR}.${PROJECT_PATCH})
  elseif(git_version MATCHES
                 "^(0|[1-9][0-9]*)[.](0|[1-9][0-9]*)")
    set(PROJECT_MAJOR "${CMAKE_MATCH_1}")
    set(PROJECT_MINOR "${CMAKE_MATCH_2}")
    set(PROJECT_VERSION ${PROJECT_MAJOR}.${PROJECT_MINOR})
  elseif(git_version MATCHES
                 "^(0|[1-9][0-9]*)")
    set(PROJECT_MAJOR "${CMAKE_MATCH_1}")
    set(PROJECT_VERSION ${PROJECT_MAJOR})
  endif()
endif()
