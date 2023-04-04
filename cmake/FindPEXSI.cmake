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
    PATH_SUFFIXES "include"
)

find_library(PEXSI_LIBRARY
    NAMES pexsi_linux_release_v2.0
    HINTS ${PEXSI_DIR}
    PATH_SUFFIXES "lib"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PEXSI DEFAULT_MSG PEXSI_INCLUDE_DIR PEXSI_LIBRARY)

if(PEXSI_FOUND)
    set(PEXSI_LIBRARIES ${PEXSI_LIBRARY})
    set(PEXSI_INCLUDE_DIR ${PEXSI_INCLUDE_DIR})

    if(NOT TARGET PEXSI::PEXSI)
        add_library(PEXSI::PEXSI UNKNOWN IMPORTED)
        set_target_properties(PEXSI::PEXSI PROPERTIES
            IMPORTED_LINK_INTERFACE_LANGUAGES "C"
            IMPORTED_LOCATION "${PEXSI_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${PEXSI_INCLUDE_DIR}")
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${PEXSI_INCLUDE_DIR})