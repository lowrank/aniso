#include <iostream>
#include "Aniso.h"
#include "utility/config.h"
#include "gmres.h"

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

        timer.tic("Cache");
        aniso.runKernelsCache(scaledFunction, output);
        timer.toc();

        timer.tic("Removal");
        aniso.nearRemoval(scaledFunction, output);
        timer.toc();

        timer.tic("NearAddOn Cache");
        aniso.refineAddOnCache(scaledFunction, output);
        timer.toc();

        timer.tic("SingularAddOn");
        aniso.singularAdd(unscaledCoefficient, output);
        timer.toc();

        return output;
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

        timer.tic("SingularAddOn");
        aniso.singularAdd(unscaledCoefficient, output);
        timer.toc();

        return output;
    };

    Vector charge(aniso.numberOfNodes);
    setValue(charge, 1.0);

    Vector output = cache_mapping(charge);

    for (int i = 0; i < aniso.numberOfNodes; ++i) {
        output(i) = output(i) / (2 * M_PI);
    }

    auto forwardOperator = [&](Vector& scaledCharge) {
        Vector load(aniso.numberOfNodes);
        // there are faster ways
        for (int i = 0; i < aniso.numberOfNodes; ++i) {
            load(i) = scaledCharge(i) * aniso.sigma_s(i);
        }
        Vector scatter = apply_mapping(load);

        // there are faster ways
        for (int i = 0; i < aniso.numberOfNodes; ++i) {
            scatter(i) = scaledCharge(i) - scatter(i) / (2 * M_PI);
        }

        return scatter;
    };

    Vector x(aniso.numberOfNodes);
    setValue(x, 0.);
    GMRES(forwardOperator, x, output, 20, 400, 1e-14);

    aniso.displayKernelCacheSize();

}