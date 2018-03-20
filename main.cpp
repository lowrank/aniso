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
    Vector charge(aniso.numberOfNodes);

    auto source_function = [&](scalar_t x, scalar_t y) {
        scalar_t d = SQR(x - 0.5) + SQR(y - 0.5);
        return exp(-25 * d);
    };

    auto scattering_function = [&](scalar_t x, scalar_t y) {
        return 10 * 0.5 * (1 - cos(2 * M_PI * x)) ;
    };

    auto total_function = [&](scalar_t x, scalar_t y) {
        return 10 * 0.5 * (1 - cos(2 * M_PI * x)) + 0.2;
    };

    for (int i = 0; i < aniso.nodes.size(); ++i) {
        charge(i) = source_function(aniso.nodes[i].x, aniso.nodes[i].y);
        aniso.sigma_t(i) = total_function(aniso.nodes[i].x, aniso.nodes[i].y) ;
        aniso.sigma_s(i) = scattering_function(aniso.nodes[i].x, aniso.nodes[i].y) ;
    }

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
        Vector output_non;
        vector<Vector> unscaledCoefficientCache;
        for (int i = 0; i < aniso.nodes.size(); ++i) {
            unscaledFunctionCache(i) = charge(i);
            scaledFunctionCache(i) = charge(i) * aniso.weights[i];
        }

        timer.tic("interpolate source");
        aniso.interpolation(unscaledFunctionCache, unscaledCoefficientCache);
        timer.toc();

        timer.tic("Non Singular Cache");
        aniso.runKernelsCache(scaledFunctionCache, output_non);
        timer.toc();

        timer.tic("Cache");
        aniso.runKernelsCacheSing(scaledFunctionCache, output_cache);
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

        daxpy(1.0, output_non, output_cache);
        dscal(M_1_PI/2.0, output_cache);
        return output_cache;
    };

    auto apply_mapping = [&](Vector& charge) {
        assert(charge.row() == aniso.nodes.size());
        Vector scaledFunction(aniso.numberOfNodes);
        Vector unscaledFunction(aniso.numberOfNodes);
        Vector output;
        Vector output_s;

        vector<Vector> unscaledCoefficient;
        for (int i = 0; i < aniso.nodes.size(); ++i) {
            unscaledFunction(i) = charge(i);
            scaledFunction(i) = charge(i) * aniso.weights[i];
        }

        timer.tic("interpolate source");
        aniso.interpolation(unscaledFunction, unscaledCoefficient);
        timer.toc();

        timer.tic("Non Singular");
        aniso.runKernelsFast(scaledFunction, output_s);
        timer.toc();

        timer.tic("Apply Fast");
        aniso.runKernelsFastSing(scaledFunction, output);
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

        daxpy(1.0, output_s, output);
        dscal(M_1_PI/2.0, output);

        return output;
    };



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
    GMRES(forwardOperator, x, rhs, 20, 400, 1e-12);

    if (atoi(cfg.options["IO"].c_str())) {
        write_to_csv(aniso.nodes, "points.csv", " ");
        write_to_csv(x, "result.csv");
    }

    aniso.displayKernelCacheSize();

}