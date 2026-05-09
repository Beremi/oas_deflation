#ifndef _OPENMP_UTILS_H
#define _OPENMP_UTILS_H

#include <cstddef>

#ifdef _OPENMP
 #include <omp.h>
#endif

namespace OmpUtils {

inline int maxThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

inline int threadNum() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

inline bool shouldParallelize(std :: size_t work, std :: size_t minWork = 256) {
    return maxThreads() > 1 && work >= minWork;
}

inline int usableThreads(std :: size_t work, std :: size_t minWork = 256) {
    if ( !shouldParallelize(work, minWork) ) {
        return 1;
    }
    return maxThreads();
}

} // namespace OmpUtils

#endif /* _OPENMP_UTILS_H */
