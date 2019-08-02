// Based on example code by Guillaume Cottenceau from libpng

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "../../3rdparty/libigl/external/stb/stb_image.h"

unsigned char * read_png_file(std::string file_name, int &width, int &height, int &channels) {
    unsigned char *image = stbi_load(file_name.c_str(), &width, &height, &channels, STBI_grey);

    return image;
}