
#ifndef __MARA_CONFIG_HEADER__
#define __MARA_CONFIG_HEADER__

#define __MARA_BASE_VERSION "1.0"
#ifdef USE_MPI
#define __MARA_USE_MPI 1
#else
#define __MARA_USE_MPI 0
#endif // USE_MPI

#ifdef __INTEL_COMPILER
#define Mara_isinf_cxx(x) isinf(x)
#define Mara_isnan_cxx(x) isnan(x)
#else
#define Mara_isinf_cxx(x) std::isinf(x)
#define Mara_isnan_cxx(x) std::isnan(x)
#endif // __INTEL_COMPILER

#endif // __MARA_CONFIG_HEADER__
