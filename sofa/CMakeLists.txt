cmake_minimum_required(VERSION 2.8.3)

file(GLOB_RECURSE  SOURCES 20170420/c/*.c)

add_library(sofa STATIC ${SOURCES})


