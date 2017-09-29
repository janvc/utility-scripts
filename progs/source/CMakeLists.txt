find_package(Boost COMPONENTS system filesystem program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Armadillo REQUIRED)
find_package(FFTW REQUIRED)

set(convolute_SOURCES
    convolute.cpp
)

set(gauss_viban_cpp_SOURCES
    gauss_viban.cpp
    GaussFchk.cpp
    utilities.cpp
    vibrationalanalysis.cpp
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
    vibrationalanalysis.cpp
	utilities.cpp
)

set(test_cubic_SOURCES
    test_cubic.cpp
    GaussFchk.cpp
    utilities.cpp
)

set(anharmgaussjob_SOURCES
    anharmgaussjob.cpp
    GaussFchk.cpp
    utilities.cpp
    vibrationalanalysis.cpp
)

set(gva2_SOURCES
    gva2.cpp
    vibrationalanalysis.cpp
    GaussFchk.cpp
    utilities.cpp
    mkduschclass.cpp
)

set(mkivrcalc_SOURCES
    mkivrcalc.cpp
    vibrationalanalysis.cpp
    GaussFchk.cpp
    utilities.cpp
)

set(mctdh_fft_SOURCES
    mctdh_fft.cpp
)

set(gau2fcc_SOURCES
    gau2fcc.cpp
    GaussFchk.cpp
    vibrationalanalysis.cpp
    utilities.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008 -fall-intrinsics -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -Og")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ggdb3 -Og")

set(Lapack_LIBRARIES lapack)

add_executable(convolute       ${convolute_SOURCES})
add_executable(gauss_viban_cpp ${gauss_viban_cpp_SOURCES})
add_executable(gauss_viban_f   ${gauss_viban_f_SOURCES})
add_executable(mctdh_spec_fft  ${mctdh_spec_fft_SOURCES})
add_executable(mctdh_specinput ${mctdh_specinput_SOURCES})
add_executable(mkdusch         ${mkdusch_SOURCES})
add_executable(test_cubic      ${test_cubic_SOURCES})
add_executable(anharmgaussjob  ${anharmgaussjob_SOURCES})
add_executable(gva2            ${gva2_SOURCES})
add_executable(mkivrcalc       ${mkivrcalc_SOURCES})
add_executable(mctdh_fft       ${mctdh_fft_SOURCES})
add_executable(gau2fcc         ${gau2fcc_SOURCES})

target_link_libraries(anharmgaussjob ${Boost_LIBRARIES})
target_link_libraries(convolute ${Boost_LIBRARIES})
target_link_libraries(gauss_viban_cpp ${Boost_LIBRARIES})
target_link_libraries(mctdh_specinput ${Boost_LIBRARIES})
target_link_libraries(mkdusch ${Boost_LIBRARIES})
target_link_libraries(test_cubic ${Boost_LIBRARIES})
target_link_libraries(gauss_viban_f ${Lapack_LIBRARIES})
target_link_libraries(gva2 ${Boost_LIBRARIES})
target_link_libraries(mkivrcalc ${Boost_LIBRARIES})
target_link_libraries(mctdh_fft ${Boost_LIBRARIES} ${FFTW_LIBRARIES})
target_link_libraries(gau2fcc ${Boost_LIBRARIES})

install(TARGETS convolute DESTINATION bin)
install(TARGETS gauss_viban_cpp DESTINATION bin)
install(TARGETS gauss_viban_f DESTINATION bin)
install(TARGETS mctdh_spec_fft DESTINATION bin)
install(TARGETS mctdh_specinput DESTINATION bin)
install(TARGETS mkdusch DESTINATION bin)
install(TARGETS test_cubic DESTINATION bin)
install(TARGETS gva2 DESTINATION bin)
install(TARGETS anharmgaussjob DESTINATION bin)
install(TARGETS mkivrcalc DESTINATION bin)
install(TARGETS mctdh_fft DESTINATION bin)
install(TARGETS gau2fcc DESTINATION bin)

