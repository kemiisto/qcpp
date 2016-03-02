cmake_minimum_required(VERSION 2.8)
project(01)
enable_language(Fortran)

if(NOT MINGW)
  find_package(LAPACK REQUIRED)
else()
  find_package(OpenBLAS)
  if(OpenBLAS_FOUND)
    MESSAGE(STATUS "OpenBLAS found.")
    set(LAPACK_LIBRARIES "libopenblas.dll.a")
  else()
    set(LAPACK_DIR "c:/developer/libraries/lapack/3.5.0/mingw-builds-5.2.0-64/debug/lib")
    set(LAPACK_LIBRARIES "${LAPACK_DIR}/liblapack.a;${LAPACK_DIR}/libblas.a")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static")
  endif()
endif()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(DEBUG_FLAGS -g -O0)
  set(RELEASE_FLAGS -O2)
  set(STANDARD_FLAGS -std=f2008 -pedantic-errors)
  # -Wall currently includes (as of GCC 5.2.0)
  # -Waliasing
  # -Wampersand
  # -Wc-binding-type
  # -Wconversion
  # -Winteger-division
  # -Wintrinsic-shadow
  # -Wintrinsics-std
  # -Wline-truncation
  # -Wno-tabs
  # -Wreal-q-constant
  # -Wsurprising
  # -Wtarget-lifetime
  # -Wunused

  # -Wextra currently includes (as of GCC 5.2.0)
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
    set(DEBUG_FLAGS -g -O0)
    set(RELEASE_FLAGS -O2)
    set(STANDARD_FLAGS -stand f03)
    set(WARNING_FLAGS -warn all -warn errors)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
    set(DEBUG_FLAGS -g -O0 -Mbounds)
    set(RELEASE_FLAGS -fast -Mipa=fast,inline)
    #set(STANDARD_FLAGS -stand f03)
    set(WARNING_FLAGS -Minform=inform)
endif()

set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "")
set(CMAKE_Fortran_FLAGS_RELEASE "")

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
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG: ${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")

add_executable(project_01 ../fcl/fcl_constants.f90 ../fcl/fcl_kinds.f90 ../fcl/fcl_lapack.f90 ../fcl/fcl_vecmath_vector_3d.F90 chem.f90 project_01.f90)
target_link_libraries(project_01 ${LAPACK_LIBRARIES})

add_executable(project_02 ../fcl/fcl_constants.f90 ../fcl/fcl_kinds.f90 ../fcl/fcl_lapack.f90 ../fcl/fcl_vecmath_vector_3d.F90 chem.f90 project_02.f90)
target_link_libraries(project_02 ${LAPACK_LIBRARIES})

add_executable(project_03 ../fcl/fcl_kinds.f90 project_03.f90)

