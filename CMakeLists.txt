cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
project(1m_TCS)

add_subdirectory(SFML-2.4.1)
include_directories(SFML-2.4.1/include)

add_subdirectory(Serial)
include_directories(Serial/include)


include_directories(eigen)

add_subdirectory(sofa)
include_directories(sofa/20170420/c/src)

add_executable(1m_TCS src/main.cpp src/pointing.cpp src/pointing.h src/site.h )

target_link_libraries(1m_TCS sfml-system sfml-window sfml-graphics serial sofa)