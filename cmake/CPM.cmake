set(CPM_VERSION 0.40.2)
set(CPM_DOWNLOAD_URL https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_VERSION}/CPM.cmake)
set(CPM_DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_VERSION}.cmake)

if(NOT EXISTS ${CPM_DOWNLOAD_DIR})
  message(STATUS "Downloading CPM.cmake to ${CPM_DOWNLOAD_DIR}")
  file(DOWNLOAD ${CPM_DOWNLOAD_URL} ${CPM_DOWNLOAD_DIR})
endif()

include(${CPM_DOWNLOAD_DIR})
