if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall
  $<$<COMPILE_LANG_AND_ID:C,AppleClang,Clang>:-Wno-unused-but-set-variable>)
endif()

# disable nuisance warnings from Visual Studio
add_compile_definitions($<$<BOOL:${MSVC}>:_CRT_SECURE_NO_WARNINGS>)
