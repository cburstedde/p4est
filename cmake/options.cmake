option(enable_p6est "build p6est" on)
option(enable_p8est "build p8est" on)

option(P4EST_BUILD_TESTING "build p4est testing" on)

option(enable-file-deprecated "use deprecated data file format" off)

option(vtk_binary "VTK binary interface" on)
if(vtk_binary)
  set(P4EST_ENABLE_VTK_BINARY 1)
endif()

# Necessary for shared library with Visual Studio / Windows oneAPI
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS true)

# --- auto-ignore build directory
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  file(GENERATE OUTPUT .gitignore CONTENT "*")
endif()
