#include "rf_pipelines_internals.hpp"

#ifdef HAVE_PNG
#include <png.h>
#endif

using namespace std;

namespace rf_pipelines {
#if 0
};  // pacify emacs c-mode!
#endif

#ifndef HAVE_PNG

void write_rgb8_png(const string &filename, uint8_t *rgb, int m, int n, bool ymajor, bool ytop_to_bottom)
{
    throw runtime_error("write_rgb8_png() was called, but this version of librf_pipelines wasn't compiled with libpng");
}

#else // HAVE_PNG


// -------------------------------------------------------------------------------------------------
//
// png_writer: RAII wrapper for libpng.


class png_writer {
public:
    // @rgb should be an array of shape (ny,nx,3) where the last index is "rgb".
    // The 'ytop_to_bottom' flag controls whether the outer index ("y") runs from
    // top to bottom in the image, or from bottom to top.

    void write_rgb8(const string &filename, uint8_t *rgb_ymajor, int nx, int ny, bool ytop_to_bottom);

    png_writer() { }
    ~png_writer();

    // Moveable but noncopyable
    png_writer(png_writer &&src);
    png_writer(const png_writer &) = delete;
    png_writer& operator=(png_writer &&src);
    png_writer& operator=(const png_writer &) = delete;

protected:
    void move(png_writer &src);
    void deallocate();

    png_structp png_ptr = nullptr;
    png_infop info_ptr = nullptr;
    FILE *fp = nullptr;
};


png_writer::~png_writer()
{
    this->deallocate();
}

png_writer::png_writer(png_writer &&src)
{
    this->move(src);
}

png_writer& png_writer::operator=(png_writer &&src)
{
    if (this != &src)
	this->move(src);
    return *this;
}


void png_writer::deallocate()
{
    // close file first
    if (fp) {
	fclose(fp);
	fp = NULL;
    }
    
    // The libpng documentation doesn't really explain how to clean up its state.
    //
    // I dug into the source code and decided that one call to png_destroy_write_struct()
    // destroys both the png_ptr and the info_ptr, and does the right thing if only the
    // png_ptr is non-NULL.

    if (png_ptr || info_ptr) {
	png_destroy_write_struct(&png_ptr, &info_ptr);
	png_ptr = NULL;
	info_ptr = NULL;
    }
}


void png_writer::move(png_writer &src)
{
    this->deallocate();

    this->png_ptr = src.png_ptr;
    this->info_ptr = src.info_ptr;
    this->fp = src.fp;

    src.png_ptr = nullptr;
    src.info_ptr = nullptr;
    src.fp = nullptr;
}


void png_writer::write_rgb8(const string &filename, uint8_t *rgb_ymajor, int nx, int ny, bool ytop_to_bottom)
{
    // Just in case any state remains from a previous call
    // FIXME (low-priority, just curious): is it necessary to reallocate the 'png_structp' and 'png_infop' between writes?
    this->deallocate();

    this->png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    rf_assert(png_ptr);

    this->info_ptr = png_create_info_struct(png_ptr);
    rf_assert(info_ptr);

    fp = fopen(filename.c_str(), "wb");

    if (!fp) {
        cerr << "rf_pipelines: couldn't open file " << filename << endl;
	this->deallocate();
        return;
    }
    
    if (setjmp(png_jmpbuf(png_ptr))) {
	cerr << "rf_pipelines: libpng internal failure (setjmp(png_jmpbuf(...)))\n";
	this->deallocate();
	return;
    }

    png_init_io(png_ptr, fp);

    png_set_IHDR(png_ptr, info_ptr, nx, ny,
		 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    uint8_t *rgb_base = ytop_to_bottom ? &rgb_ymajor[0] : &rgb_ymajor[(ny-1)*(3*nx)];
    int rgb_ystride = ytop_to_bottom ? (3*nx) : (-3*nx);

    for (int i = 0; i < ny; i++)
	png_write_row(png_ptr, reinterpret_cast<png_byte *> (rgb_base + i * rgb_ystride));

    png_write_end(png_ptr, NULL);

    this->deallocate();
    cout << "rf_pipelines: wrote " << filename << endl;
}


// -------------------------------------------------------------------------------------------------


// The 'rgb' array should have shape (m,n,3).
//
// The 'ymajor' flag controls whether the major (length-m) index of the array is the y-axis of the
// plots, i.e. the plot dimensions are either (n,m) or (m,n) depending on whether ymajor is true or false.

void write_rgb8_png(const string &filename, uint8_t *rgb, int m, int n, bool ymajor, bool ytop_to_bottom)
{
    // temporary buffer, if transpose is needed.
    vector<uint8_t> tmp;

    if (!ymajor) {
	tmp.resize(n*m*3);

	for (int i = 0; i < m; i++) {
	    for (int j = 0; j < n; j++) {
		tmp[j*m*3 + i*3] = rgb[i*n*3 + j*3];
		tmp[j*m*3 + i*3 + 1] = rgb[i*n*3 + j*3 + 1];
		tmp[j*m*3 + i*3 + 2] = rgb[i*n*3 + j*3 + 2];
	    }
	}
	
	rgb = &tmp[0];
	std::swap(m, n);
    }

    png_writer pw;
    pw.write_rgb8(filename, rgb, n, m, ytop_to_bottom);
}

#endif  // HAVE_PNG


}   // namespace rf_pipelines
