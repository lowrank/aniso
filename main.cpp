#include <iostream>
#include "Aniso.h"
int main() {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif
    Aniso aniso(128, 1, 1, 0.8, 64);

    std::cout << aniso.anisotropy << std::endl;

    // load function.
    setValue(aniso.sigma_t, 40.2);
    setValue(aniso.sigma_s, 40.0);

    aniso.interpolation();

    Vector f(aniso.numberOfNodes);
    setValue(f, 1.);
    aniso.makeKernels();
    aniso.runKernels(f);
}