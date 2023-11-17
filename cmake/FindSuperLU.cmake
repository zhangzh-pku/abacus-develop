###############################################################################
# - Find cereal
# Find the native cereal headers.
#
#  CEREAL_FOUND - True if cereal is found.
#  CEREAL_INCLUDE_DIR - Where to find cereal headers.

# find_path(SuperLU_INCLUDE_DIR
#     NAMES *.h
#     HINTS ${SuperLU_DIR}
#     PATH_SUFFIXES "include"
# )

find_library(SuperLU_LIBRARY
    NAMES libsuperlu_dist.a
    HINTS ${SuperLU_DIR}
    PATH_SUFFIXES "lib"
)

# Handle the QUIET and REQUIRED arguments and
# set Cereal_FOUND to TRUE if all variables are non-zero.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SuperLU DEFAULT_MSG SuperLU_LIBRARY)

# Copy the results to the output variables and target.
mark_as_advanced(SuperLU_LIBRARY)

