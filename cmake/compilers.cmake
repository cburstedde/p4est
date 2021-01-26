if(CMAKE_C_COMPILER_ID STREQUAL GNU)
  add_compile_options(-fmax-errors=3)  # avoid extremely long error printouts
endif()
