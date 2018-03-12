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

    Profiler timer;

    vector<scalar_t > a,b,c;

    aniso.duffy_transform({0, 0,  1. , 0  , 1. ,   1.}, a,b,c);

    for (auto cc : c) {
        std::cout << cc << std::endl;
    }

    // load function.
    setValue(aniso.sigma_t, 40.2);
    setValue(aniso.sigma_s, 40.0);

    aniso.interpolation();

    Vector f(aniso.numberOfNodes);
    Vector h(aniso.numberOfNodes);
    vector<Vector> h_coeff;

    setValue(f, 1.);

    for (int i = 0; i < aniso.nodes.size(); ++i) {
        h(i) = f(i);
        f(i) = f(i) * aniso.weights[i];
    }

    aniso.interpolation(h, h_coeff);
    std::cout << h_coeff[0] << std::endl;

    timer.tic("?");
    aniso.singPrecompute();
    timer.toc();

    scalar_t s = 0.;
    for (int i = 0; i < aniso.singW[0].size(); ++i) {
        s += aniso.singW[0][i];
    }

    std::cout << s << std::endl;

    aniso.makeKernels();



//    timer.tic("Normal");
//    aniso.runKernels(f);
//    timer.toc();

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

    timer.tic("NearAddOn Cache");
    aniso.refineAddOnCache(f);
    timer.toc();

    timer.tic("NearAddOn Fast");
    aniso.refineAddOnFast(f);
    timer.toc();

    timer.tic("SingularAddOn");
    aniso.singularAdd(f);
    timer.toc();

    aniso.displayKernelCacheSize();

}