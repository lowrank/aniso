#include <iostream>
#include "Aniso.h"
#include "utility/config.h"

int main(int argc, char* argv[]) {
#ifdef RUN_OMP
    omp_set_num_threads(omp_get_max_threads());
#endif

    config cfg(argv[1]); cfg.print();

    auto domain_size     = atoi(cfg.options["domainSize"].c_str());
    auto quadrature_rule = atoi(cfg.options["quadRule"].c_str());
    auto kernel_size     = atoi(cfg.options["kernelSize"].c_str());
    auto anisotropy      = atof(cfg.options["g"].c_str());
    auto singular_rule   = atoi(cfg.options["singRule"].c_str());
    auto fmm_np          = atoi(cfg.options["np"].c_str());
    auto fmm_maxLevel    = atoi(cfg.options["maxLevel"].c_str());

    Aniso aniso(domain_size, quadrature_rule,  kernel_size, anisotropy, singular_rule, fmm_np, fmm_maxLevel);

    // load function.
    setValue(aniso.sigma_t, 40.2);
    setValue(aniso.sigma_s, 40.0);

    aniso.interpolation();

    Vector f(aniso.numberOfNodes);
    setValue(f, 1.);
    for (int i = 0; i < aniso.nodes.size(); ++i) {
        f(i) = f(i) * aniso.weights[i];
    }

    aniso.makeKernels();

    Profiler timer;

    std::cout << aniso.nodes[0].x << " " << aniso.nodes[0].y << std::endl;

    timer.tic("Normal");
    aniso.runKernels(f);
    timer.toc();

#ifdef BBFMM_CACHE
    timer.tic("Cache");
    aniso.runKernelsCache(f);
    timer.toc();

    timer.tic("Apply");
    aniso.runKernelsFast(f);
    timer.toc();
#endif

    timer.tic("Removal");
    aniso.nearRemoval(f);
    timer.toc();

    aniso.displayKernelCacheSize();

}