cmake_minimum_required(VERSION 2.8)
project(01)
enable_language(Fortran)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(DEBUG_FLAGS -O0)
  set(RELEASE_FLAGS -O2)
  set(STANDARD_FLAGS -std=f2008 -pedantic-errors)
  # -Wall currently includes (as of GCC 4.9.2)
  # -Waliasing
  # -Wampersand
  # -Wc-binding-type
  # -Wconversion
  # -Wintrinsic-shadow
  # -Wintrinsics-std
  # -Wline-truncation
  # -Wno-tabs
  # -Wreal-q-constant
  # -Wsurprising
  # -Wtarget-lifetime
  # -Wunused

  # -Wextra currently includes (as of GCC 4.9.2)
  # -Wcompare-reals
  # -Wunused-parameter
  set(WARNING_FLAGS -Wall -Wextra -Werror
    # -Warray-temporaries
    -Wcharacter-truncation
    -Wconversion-extra
    -Wimplicit-interface
    -Wimplicit-procedure
    -Wunderflow
    -fbounds-check
    -fbacktrace
    -ffpe-trap=zero
    -fimplicit-none
    -fmax-errors=1)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set(DEBUG_FLAGS -O0)
    set(RELEASE_FLAGS -O2)
    set(STANDARD_FLAGS -stand f03)
    set(WARNING_FLAGS -warn all -warn errors)
endif()

foreach(FLAG ${STANDARD_FLAGS} ${WARNING_FLAGS})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLAG}")
endforeach()

foreach(FLAG ${DEBUG_FLAGS})
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${FLAG}")
endforeach()

foreach(FLAG ${RELEASE_FLAGS})
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${FLAG}")
endforeach()

message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")

add_executable(project_01 kinds.f90 vecmath.f90 chem.f90 project_01.f90)