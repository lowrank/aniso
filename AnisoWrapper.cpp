//
// Created by lurker on 4/7/18.
//

#include "AnisoWrapper.h"

template class Session<Aniso>;

namespace {
    MEX_DEFINE(new) (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 7);
        OutputArguments output(nlhs, plhs, 1);

        auto domain_size     =  M_Cast<double>(C_CAST(prhs[0]));
        auto quadrature_rule =  M_Cast<double>(C_CAST(prhs[1]));
        auto kernel_size     =  M_Cast<double>(C_CAST(prhs[2]));
        auto anisotropy      =  M_Cast<double>(C_CAST(prhs[3]));
        auto singular_rule   =  M_Cast<double>(C_CAST(prhs[4]));
        auto fmm_np          =  M_Cast<double>(C_CAST(prhs[5]));
        auto fmm_maxLevel    =  M_Cast<double>(C_CAST(prhs[6]));

        output.set(0, Session<Aniso>::create(new Aniso(int(*domain_size), int(*quadrature_rule),
                                                       int(*kernel_size), *anisotropy, int(*singular_rule),
                                                       int(*fmm_np), int(*fmm_maxLevel))));
    }

    MEX_DEFINE(delete)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 0);
        Session<Aniso>::destroy(input.get(0));
    }

    MEX_DEFINE(getNodes)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 1);
        auto as = Session<Aniso>::get(input.get(0));
        plhs[0] = mxCreateNumericMatrix((int32_t)(as->nodes.size()),2, mxDOUBLE_CLASS, mxREAL);
        auto n_ptr = mxGetPr(plhs[0]);
        for (int i = 0; i < as->nodes.size(); ++i) {
            n_ptr[i                   ] = as->nodes[i].x;
            n_ptr[i + as->nodes.size()] = as->nodes[i].y;
        }

    }

    MEX_DEFINE(setCoeff)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 3);
        OutputArguments output(nlhs, plhs, 0);
        auto as = Session<Aniso>::get(input.get(0));

        // get information
        auto sigma_s_ptr = M_Cast<double>(C_CAST (prhs[1]));
        auto sigma_t_ptr = M_Cast<double>(C_CAST (prhs[2]));

        assert(mxGetM(prhs[1]) == as->nodes.size());
        assert(mxGetM(prhs[2]) == as->nodes.size());

        for (int i = 0; i < as->nodes.size(); ++i) {
            as->sigma_t(i) = sigma_t_ptr[i];
            as->sigma_s(i) = sigma_s_ptr[i];
        }

        as->interpolation();

        as->singPrecompute();

        as->makeKernels();

    }


    MEX_DEFINE(cache)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 0);
        auto as = Session<Aniso>::get(input.get(0));

        as->runKernelsCache(0);

        as->runKernelsCacheSing(0);

        as->refineAddOnCache(0);

        as->singularAddCache(0);

    }

    MEX_DEFINE(mapping)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        InputArguments input(nrhs, prhs, 2);
        OutputArguments output(nlhs, plhs, 1);
        auto as = Session<Aniso>::get(input.get(0));

        auto charge = M_Cast<double>(C_CAST (prhs[1]));

        Vector scaledFunction(as->numberOfNodes);
        Vector unscaledFunction(as->numberOfNodes);
        Vector output_t;
        Vector output_s;

        plhs[0] = mxCreateNumericMatrix((int32_t)(as->nodes.size()),1, mxDOUBLE_CLASS, mxREAL);
        auto n_ptr = mxGetPr(plhs[0]);


        vector<Vector> unscaledCoefficient;
        for (int i = 0; i < as->nodes.size(); ++i) {
            unscaledFunction(i) = charge[i];
            scaledFunction(i) = charge[i] * as->weights[i];
        }

        as->interpolation(unscaledFunction, unscaledCoefficient);

        as->runKernelsFast(0,scaledFunction, output_s);

        as->runKernelsFastSing(0,scaledFunction, output_t);

        as->nearRemoval(0,scaledFunction, output_t);

        as->refineAddOnFast(0,scaledFunction, output_t);

        as->singularAddFast(0,unscaledCoefficient, output_t);

        daxpy(1.0, output_s, output_t);
        dscal(M_1_PI/2.0, output_t);


        for (int i = 0; i < as->numberOfNodes; ++i) {
            n_ptr[i] = output_t(i);
        }
    }
}
MEX_DISPATCH