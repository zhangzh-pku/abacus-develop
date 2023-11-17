###############################################################################
# - Find cereal
# Find the native cereal headers.
#
#  CEREAL_FOUND - True if cereal is found.
#  CEREAL_INCLUDE_DIR - Where to find cereal headers.

find_path(ParMETIS_INCLUDE_DIR
    NAMES metis.h parmetis.h
    HINTS ${ParMETIS_DIR}
    PATH_SUFFIXES "include"
)

find_library(METIS_LIBRARY
    NAMES metis
    HINTS ${ParMETIS_DIR}
    PATH_SUFFIXES "lib"
)

find_library(ParMETIS_LIBRARY
    NAMES parmetis
    HINTS ${ParMETIS_DIR}
    PATH_SUFFIXES "lib"
)

# print libs
# message(STATUS "ParMETIS_INCLUDE_DIR: ${ParMETIS_INCLUDE_DIR}")
# message(STATUS "ParMETIS_LIBRARY: ${ParMETIS_LIBRARY}")

# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS DEFAULT_MSG ParMETIS_LIBRARY METIS_LIBRARY)

# Copy the results to the output variables and target.
mark_as_advanced(ParMETIS_LIBRARY)

