
#ifndef __MARA_CONFIG_HEADER__
#define __MARA_CONFIG_HEADER__

#define __MARA_BASE_VERSION "1.0"

#define __MARA_INSTALL_DIR "/Users/jzrake/Work/Mara"
#define __MARA_USE_MPI 1
#define __MARA_USE_HDF5 0
#define __MARA_USE_HDF5_PAR 0
#define __MARA_USE_FFTW 0
#define __MARA_USE_GLFW 0

#ifdef __INTEL_COMPILER
#define Mara_isinf_cxx(x) isinf(x)
#define Mara_isnan_cxx(x) isnan(x)
#else
#define Mara_isinf_cxx(x) std::isinf(x)
#define Mara_isnan_cxx(x) std::isnan(x)
#endif // __INTEL_COMPILER

#endif // __MARA_CONFIG_HEADER__
