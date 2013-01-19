
#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MaraImage_HEADER__
#define __MaraImage_HEADER__

void Mara_image_write_ppm(const char *fname, const double *data, int cmap,
			  int Nx, int Ny, const double *range);
const float *Mara_image_get_colormap(int cmap);

#endif // __MaraImage_HEADER__

#ifdef __cplusplus
}
#endif
