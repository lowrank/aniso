//
// Created by lurker on 2/4/18.
//

#ifndef ANISO2_KERNELFACTORY_H
#define ANISO2_KERNELFACTORY_H

#include "bbfmm/bbfmm.h"
#include "Quadrature.h"
#include "Geometry.h"

using namespace std::placeholders;

class KernelFactory : public  Geometry{
public:
    KernelFactory(){}
    KernelFactory(int geometry_size, int geometry_degree,
                  int kernel_size, double kernel_g, int kernel_degree);
    virtual ~KernelFactory();

    // amazing, there is only real parts taking effect.
    vector<kernel> realParts;
    int kernelSize;

    Vector anisotropy;

    Quadrature lineQuadratureRule;
    Quadrature singQuadratureRule;

    Vector sigma_t;
    Vector sigma_s;

    index_t np;
    index_t maxLevel;

    // a redundant storage. O(KernelSize * numberOfNodes * 9 * RefineQuadrature = K * 4^L * N^2 * Q^2 * 9)
    vector<vector<vector<vector<vector<scalar_t > > > > > nearInteractions;

    scalar_t lineIntegral(point& p, point &q);

    void interpolation();
    void makeKernels();
    void runKernels(Vector& f);
    void runKernelsCache(Vector& f);
    void runKernelsFast(Vector& f);

    void nearRemoval(Vector& f);
    void refineAddOnCache(Vector& f);
    void refineAddOnFast(Vector& f);
    void singularAdd(Vector& f);

    vector<Vector> sigma_t_coeff;
    vector<Vector> sigma_s_coeff;


    scalar_t lineIntegral(double x0, double y0, double x1, double y1);
    scalar_t integral_helper(double x0, double y0, double x1, double y1);

    int getRow(double y) ;

    int getCol(double x) ;


};


#endif //ANISO2_KERNELFACTORY_H
