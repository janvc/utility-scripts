find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(convolute_SOURCES
    convolute.cpp
)

set(gauss_viban_cpp_SOURCES
    gauss_viban.cpp
    GaussFchk.cpp
    utilities.cpp
)

set(gauss_viban_f_SOURCES
    gauss_viban.f90
    routines.f90
)

set(mctdh_spec_fft_SOURCES
    mctdh_spec_fft.f90
)

set(mctdh_specinput_SOURCES
    mctdh_specinput.cpp
    utilities.cpp
)

set(mkdusch_SOURCES
	mkdusch.cpp
	GaussFchk.cpp
	utilities.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fall-intrinsics -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ggdb3")

set(Lapack_LIBRARIES lapack)

add_executable(convolute       ${convolute_SOURCES})
add_executable(gauss_viban_cpp ${gauss_viban_cpp_SOURCES})
add_executable(gauss_viban_f   ${gauss_viban_f_SOURCES})
add_executable(mctdh_spec_fft  ${mctdh_spec_fft_SOURCES})
add_executable(mctdh_specinput ${mctdh_specinput_SOURCES})
add_executable(mkdusch         ${mkdusch_SOURCES})

target_link_libraries(convolute ${Boost_LIBRARIES})
target_link_libraries(gauss_viban_cpp ${Boost_LIBRARIES})
target_link_libraries(mctdh_specinput ${Boost_LIBRARIES})
target_link_libraries(mkdusch ${Boost_LIBRARIES})
target_link_libraries(gauss_viban_f ${Lapack_LIBRARIES})

install(TARGETS convolute DESTINATION bin)
install(TARGETS gauss_viban_cpp DESTINATION bin)
install(TARGETS gauss_viban_f DESTINATION bin)
install(TARGETS mctdh_spec_fft DESTINATION bin)
install(TARGETS mctdh_specinput DESTINATION bin)
install(TARGETS mkdusch DESTINATION bin)
