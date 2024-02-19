# else it's an offline archive
if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/.git)
  find_package(Git REQUIRED)
endif()

function(git_submodule submod_dir)
# get/update Git submodule directory to CMake, assuming the
# Git submodule directory is a CMake project.

# EXISTS because Git submodules have .git as a file, not directory
if(NOT EXISTS ${PROJECT_SOURCE_DIR}/.git)
  message(DEBUG "${PROJECT_SOURCE_DIR} is not a Git repository, skipping submodule ${submod_dir}")
  return()
endif()

if(EXISTS ${submod_dir}/CMakeLists.txt)
  return()
endif()

execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- ${submod_dir}
WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
RESULT_VARIABLE err)

if(NOT err EQUAL 0)
  message(FATAL_ERROR "${submod_dir} Git submodule failed to retrieve.")
endif()

endfunction()
