################################################################################
# Source files
################################################################################

set(SOURCES
        hexlib.h
)

micro_add_library(hex_inflator INTERFACE)
micro_copy_headers(micro_hex_inflator ${SOURCES})
