###############################################################################
# - Find PEXSI
# Find the native PEXSI headers and libraries.
#
#  PEXSI_FOUND        - True if libpexsi is found.
#  PEXSI_LIBRARIES    - List of libraries when using libpexsi
#  PEXSI_INCLUDE_DIR - Where to find PEXSI headers.
#

find_path(PEXSI_INCLUDE_DIR
    c_pexsi_interface.h
    HINTS ${PEXSI_DIR}
    PATH_SUFFIXES "include" "include/c_pexsi_interface.h"
)

find_path(PEXSI_LIBRARY
    NAMES libpexsi_linux_release_v2.0
    HINTS ${PEXSI_DIR}
    PATH_SUFFIXES "lib"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PEXSI DEFAULT_MSG PEXSI_LIBRARY PEXSI_INCLUDE_DIR)

if(PEXSI_FOUND)
    set(PEXSI_LIBRARIES ${PEXSI_LIBRARY})
    set(PEXSI_INCLUDE_DIR ${PEXSI_INCLUDE_DIR})
endif()