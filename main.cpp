#include <iostream>
#include "Aniso.h"
#include "utility/config.h"
#include "gmres.h"
#include "matlab_io.h"

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

    // load function.
    setValue(aniso.sigma_t, 20.2);
    setValue(aniso.sigma_s, 20.0);


    timer.tic("interpolate sigma");
    aniso.interpolation();
    timer.toc();

    timer.tic("precompute singular integral");
    aniso.singPrecompute();
    timer.toc();

    timer.tic("generate kernels");
    aniso.makeKernels();
    timer.toc();

    auto cache_mapping = [&](Vector& charge) {
        assert(charge.row() == aniso.nodes.size());
        Vector scaledFunctionCache(aniso.numberOfNodes);
        Vector unscaledFunctionCache(aniso.numberOfNodes);
        Vector output_cache;

        vector<Vector> unscaledCoefficientCache;
        for (int i = 0; i < aniso.nodes.size(); ++i) {
            unscaledFunctionCache(i) = charge(i);
            scaledFunctionCache(i) = charge(i) * aniso.weights[i];
        }

        timer.tic("interpolate source");
        aniso.interpolation(unscaledFunctionCache, unscaledCoefficientCache);
        timer.toc();

        timer.tic("Cache");
        aniso.runKernelsCache(scaledFunctionCache, output_cache);
        timer.toc();

        timer.tic("Removal");
        aniso.nearRemoval(scaledFunctionCache, output_cache);
        timer.toc();

        timer.tic("NearAddOn Cache");
        aniso.refineAddOnCache(scaledFunctionCache, output_cache);
        timer.toc();

        timer.tic("SingularAddOn Cache");
        aniso.singularAddCache(unscaledCoefficientCache, output_cache);
        timer.toc();
        dscal(M_1_PI/2.0, output_cache);
        return output_cache;
    };

    auto apply_mapping = [&](Vector& charge) {
        assert(charge.row() == aniso.nodes.size());
        Vector scaledFunction(aniso.numberOfNodes);
        Vector unscaledFunction(aniso.numberOfNodes);
        Vector output;

        vector<Vector> unscaledCoefficient;
        for (int i = 0; i < aniso.nodes.size(); ++i) {
            unscaledFunction(i) = charge(i);
            scaledFunction(i) = charge(i) * aniso.weights[i];
        }

        timer.tic("interpolate source");
        aniso.interpolation(unscaledFunction, unscaledCoefficient);
        timer.toc();

        timer.tic("Apply Fast");
        aniso.runKernelsFast(scaledFunction, output);
        timer.toc();

        timer.tic("Removal");
        aniso.nearRemoval(scaledFunction, output);
        timer.toc();

        timer.tic("NearAddOn Fast");
        aniso.refineAddOnFast(scaledFunction, output);
        timer.toc();

        timer.tic("SingularAddOn Fast");
        aniso.singularAddFast(unscaledCoefficient, output);
        timer.toc();

        dscal(M_1_PI/2.0, output);

        return output;
    };

    Vector charge(aniso.numberOfNodes);
    setValue(charge, 1.0);

    Vector rhs = cache_mapping(charge);

    auto forwardOperator = [&](Vector& scaledCharge) {
        Vector in = scaledCharge;
        Vector load(aniso.numberOfNodes);

        for (int i = 0; i < aniso.numberOfNodes; ++i) {
            load(i) = in(i) * aniso.sigma_s(i);
        }
        Vector scatter = apply_mapping(load);

        daxpy(-1.0 , scatter, in);
        return in;
    };

    Vector x(aniso.numberOfNodes);
    setValue(x, 0.);
    GMRES(forwardOperator, x, rhs, 5, 50, 1e-14);


    write_to_csv(aniso.nodes, "points.csv", " ");
    write_to_csv(x, "result.csv");

    aniso.displayKernelCacheSize();

}