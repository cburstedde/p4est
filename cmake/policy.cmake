set(CMAKE_EXPORT_COMPILE_COMMANDS on)

list(APPEND CMAKE_CTEST_ARGUMENTS --output-on-failure)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Debug or Release")
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.20)
  # ninja path resolution
  cmake_policy(SET CMP0116 NEW)
  # explicit source file extensions
  cmake_policy(SET CMP0115 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
  # smarter ExternalProject
  cmake_policy(SET CMP0114 NEW)
  # make missing imported targets fail immediately
  cmake_policy(SET CMP0111 NEW)
  # better find_program
  cmake_policy(SET CMP0109 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  # saner ALIAS target policies
  cmake_policy(SET CMP0107 NEW)
  cmake_policy(SET CMP0108 NEW)
  # fix link options
  cmake_policy(SET CMP0105 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.17)
  # transitive link properties
  cmake_policy(SET CMP0099 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.16)
  # ExtProj submodule
  cmake_policy(SET CMP0097 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.15)
  # better Python find
  cmake_policy(SET CMP0094 NEW)
endif()
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.14)
  # fix IN_LIST
  cmake_policy(SET CMP0085 NEW)
  # fix PIE
  cmake_policy(SET CMP0083 NEW)
  # fix add_subdirectory
  cmake_policy(SET CMP0082 NEW)
endif()
