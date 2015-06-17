cmake_minimum_required(VERSION 2.6)
project(masking C CXX)

set(CMAKE_BUILD_TYPE Release)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(masking SHARED src/masking.cpp src/spectral.cpp)
add_executable(masking-bin src/main.cpp src/readimage.cpp)
target_link_libraries(masking-bin masking)

#install libraries, binaries and header files
install(TARGETS masking DESTINATION lib)
install(FILES src/spectral.h src/masking.h DESTINATION include/masking/)
