#include "rf_pipelines_internals.hpp"

using namespace std;

namespace rf_pipelines {
#if 0
}  // emacs pacifier
#endif


zoomable_tileset::zoomable_tileset(const zoomable_tileset::initializer &ini_params_, const uint8_t background_rgb_[3]) :
    ini_params(ini_params_)
{
    this->background_rgb[0] = background_rgb_[0];
    this->background_rgb[1] = background_rgb_[1];
    this->background_rgb[2] = background_rgb_[2];

    if (!ini_params.is_initialized())
	throw runtime_error("zoomable_tileset: all fields { img_prefix, img_nzoom, img_nx, img_ny, downsample_nt } must be initialized");
}


bool zoomable_tileset::initializer::is_initialized() const
{
    return (img_prefix.size() > 0) && (img_nzoom > 0) && (img_nx > 0) && (img_ny > 0) && (downsample_nt > 0);
}

bool zoomable_tileset::initializer::is_uninitialized() const
{
    return (img_prefix.size() == 0) && (img_nzoom == 0) && (img_nx == 0) && (img_ny == 0) && (downsample_nt == 0);
}

bool zoomable_tileset::initializer::is_valid() const
{
    return is_initialized() || is_uninitialized();
}

void zoomable_tileset::initializer::reset()
{
    this->img_prefix = string();
    this->img_nzoom = 0;
    this->img_nx = 0;
    this->img_ny = 0;
    this->downsample_nt = 0;
}


}  // namespace rf_pipelines
