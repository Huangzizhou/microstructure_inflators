################################################################################
# Source files
################################################################################

set(SOURCES
        hexlib.h
)

add_library(hex_inflator INTERFACE)
micro_copy_headers(hex_inflator ${SOURCES})
