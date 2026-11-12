# bundle_runtime.cmake: copies the DLLs that gshtools._core needs next to the extension (Windows only).
#
# Inputs (set with -D in script mode, or with set() before include()):
#   GSH_PYD         the extension file. Default: the installed gshtools/_core*.pyd under CMAKE_INSTALL_PREFIX.
#   GSH_DEST        destination directory. Default: the directory of GSH_PYD.
#   GSH_EXTRA_DIRS  directories searched for DLLs that are not on PATH (empty when no third-party DLL is linked).
#
# Copied: the runtimes libomp*, vcomp*, msvcp*, vcruntime*
# (the OpenMP and MSVC runtime that clang-cl links; libomp140.x86_64.dll lives in the Windows directory on this
# toolchain, so it is selected by name and not by directory). Not copied: python3*.dll (provided by the interpreter),
# api-ms-*/ext-ms-* API sets and the other Windows system DLLs.

cmake_minimum_required(VERSION 3.24)
if(POLICY CMP0207)
  cmake_policy(SET CMP0207 NEW)  # file(GET_RUNTIME_DEPENDENCIES) normalises paths before matching regexes
endif()

if(NOT CMAKE_HOST_WIN32)
  return()
endif()

if(NOT DEFINED GSH_PYD)
  file(GLOB _gsh_pyd_candidates "${CMAKE_INSTALL_PREFIX}/gshtools/_core*.pyd")
  list(GET _gsh_pyd_candidates 0 GSH_PYD)
endif()
if(NOT DEFINED GSH_DEST)
  get_filename_component(GSH_DEST "${GSH_PYD}" DIRECTORY)
endif()

file(GET_RUNTIME_DEPENDENCIES
  LIBRARIES "${GSH_PYD}"
  DIRECTORIES ${GSH_EXTRA_DIRS}
  PRE_EXCLUDE_REGEXES "^api-ms-" "^ext-ms-" "^python3"
  RESOLVED_DEPENDENCIES_VAR _gsh_resolved
  UNRESOLVED_DEPENDENCIES_VAR _gsh_unresolved)

foreach(_dep IN LISTS _gsh_resolved)
  get_filename_component(_dep_name "${_dep}" NAME)
  string(TOLOWER "${_dep_name}" _dep_lower)
  if(_dep_lower MATCHES "^(libomp|vcomp|msvcp|vcruntime)")
    file(COPY "${_dep}" DESTINATION "${GSH_DEST}")
    message(STATUS "gshtools: bundled runtime ${_dep_name}")
  endif()
endforeach()

# The OpenMP runtime (libomp140.x86_64.dll) is not resolved by file(GET_RUNTIME_DEPENDENCIES) on this toolchain
# although it lives in the system directory, so it is copied from there by name.
set(_gsh_sys_dir "$ENV{SystemRoot}/System32")
foreach(_dep IN LISTS _gsh_unresolved)
  string(TOLOWER "${_dep}" _dep_lower)
  if(_dep_lower MATCHES "^(libomp|vcomp)")
    # Search the system directory, then the Visual Studio redistributable folders (the OpenMP runtime of clang-cl is
    # shipped in Microsoft.VC*.OpenMP.LLVM, under debug_nonredist on current Visual Studio; a CI runner may not have it in System32).
    file(GLOB _gsh_omp_candidates
      "${_gsh_sys_dir}/${_dep}"
      "$ENV{VCToolsRedistDir}/x64/Microsoft.VC*.OpenMP.LLVM/${_dep}"
      "$ENV{VCToolsRedistDir}/debug_nonredist/x64/Microsoft.VC*.OpenMP.LLVM/${_dep}"
      "$ENV{VCINSTALLDIR}/Redist/MSVC/*/x64/Microsoft.VC*.OpenMP.LLVM/${_dep}"
      "$ENV{VCINSTALLDIR}/Redist/MSVC/*/debug_nonredist/x64/Microsoft.VC*.OpenMP.LLVM/${_dep}")
    if(_gsh_omp_candidates)
      list(GET _gsh_omp_candidates 0 _gsh_omp_found)
      file(COPY "${_gsh_omp_found}" DESTINATION "${GSH_DEST}")
      message(STATUS "gshtools: bundled runtime ${_dep} (from ${_gsh_omp_found})")
    else()
      message(WARNING "gshtools: runtime ${_dep} not found in the system or Visual Studio redistributable folders; the extension will need it installed")
    endif()
  endif()
endforeach()
