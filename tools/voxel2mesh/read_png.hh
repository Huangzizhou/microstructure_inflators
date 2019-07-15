// Based on example code by Guillaume Cottenceau from libpng

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include <png.h>

png_bytep * read_png_file(std::string file_name, size_t &width, size_t &height)
{
    char header[8];    // 8 is the maximum size that can be checked

    FILE *fp = fopen(file_name.c_str(), "rb");
    if (!fp)
        throw std::runtime_error("PNG file could not be opened for reading");
    fread(header, 1, 8, fp);
    if (png_sig_cmp((png_const_bytep) header, 0, 8))
        throw std::runtime_error("File could not be recognized as a PNG file");

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        throw std::runtime_error("png_create_read_struct failed");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        throw std::runtime_error("png_create_info_struct failed");

    if (setjmp(png_jmpbuf(png_ptr)))
        throw std::runtime_error("Error during init_io");

    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    png_byte color_type = png_get_color_type(png_ptr, info_ptr);
    png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);

    png_read_update_info(png_ptr, info_ptr);

    /* read file */
    if (setjmp(png_jmpbuf(png_ptr)))
        throw std::runtime_error("Error during read_image");

    png_bytep * row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    for (size_t y=0; y<height; y++)
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));

    png_read_image(png_ptr, row_pointers);

    fclose(fp);

    //std::cout << "PNG code: " << (int) png_get_color_type(png_ptr, info_ptr) << std::endl;
    //std::cout << "Bit depth: " << (int) bit_depth << std::endl;

    if (png_get_color_type(png_ptr, info_ptr) != PNG_COLOR_TYPE_GRAY)
        throw std::runtime_error("color_type of input file must be PNG_COLOR_TYPE_GRAY");

    if (bit_depth != 1)
        throw std::runtime_error("bit depth should be 1");

    return row_pointers;
}