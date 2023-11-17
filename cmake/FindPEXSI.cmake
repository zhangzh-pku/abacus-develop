###############################################################################
# - Find cereal
# Find the native cereal headers.
#
#  CEREAL_FOUND - True if cereal is found.
#  CEREAL_INCLUDE_DIR - Where to find cereal headers.

find_path(PEXSI_INCLUDE_DIR
    NAMES c_pexsi_interface.h
    HINTS ${PEXSI_DIR}
    PATH_SUFFIXES "include"
)

find_library(PEXSI_LIBRARY
    NAMES pexsi
    HINTS ${PEXSI_DIR}
    PATH_SUFFIXES "lib"
)

# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PEXSI DEFAULT_MSG PEXSI_LIBRARY PEXSI_INCLUDE_DIR)

# Copy the results to the output variables and target.
mark_as_advanced(PEXSI_LIBRARY PEXSI_INCLUDE_DIR)

