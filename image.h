
#ifdef __cplusplus
extern "C" {
#endif

#ifndef __MaraImage_HEADER__
#define __MaraImage_HEADER__

void Mara_write_ppm(const char *fname, const double *data, int cmap,
		    int Nx, int Ny, const double *range);

#endif // __MaraImage_HEADER__

#ifdef __cplusplus
}
#endif
