include(CheckCCompilerFlag)

# --- compiler options

check_c_compiler_flag(-Wall _has_wall)
if(_has_wall)
  add_compile_options(-Wall)
else()
  check_c_compiler_flag(/Wall _has_msvc_wall)
  if(_has_msvc_wall)
    add_compile_options(/Wall)
  endif()
endif()


# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()
