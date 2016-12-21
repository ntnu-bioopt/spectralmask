cmake_minimum_required(VERSION 2.6)
project(masking C CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR})

add_library(masking SHARED src/masking.cpp src/spectral.cpp)
add_executable(masking-bin src/main.cpp src/readimage.cpp)
target_link_libraries(masking-bin masking)

set(REFLECTANCE_MASKING_SPECTRA_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/reflectance_spectra/")
set(TRANSMITTANCE_MASKING_SPECTRA_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/transmittance_spectra/")
configure_file(src/masking.h.in masking.h @ONLY)

install(TARGETS masking DESTINATION lib)
install(FILES src/spectral.h masking.h DESTINATION include/masking/)
