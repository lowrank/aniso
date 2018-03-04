//
// Created by lurker on 2/4/18.
//

/*
 * Geometry creates object to store the domain as unit square.
 *
 * 1. The unit square is divided into N^2 sub-squares.
 *
 * 2. On each sub-square, a uniform specific quadrature rule is applied.
 *
 * 3. The resulting points are stored in the output as std::vector<point>.
 *
 * 4.
 */

#ifndef ANISO2_GEOMETRY_H
#define ANISO2_GEOMETRY_H

#include "bbfmm/utils.h"
#include "bbfmm/linalg.h"
#include "bbfmm/bbfmm.h"
#include "Quadrature.h"
#include "utility/Profiler.h"

using namespace bbfmm;

class Geometry  {
public:
    Geometry(){}
    Geometry(int size, int degree);
    virtual ~Geometry();

    void makeLegendreMatrix(Matrix  &K, int N, vector<scalar_t> &x, vector<scalar_t> &y, vector<scalar_t> &w);

    //@public members
    vector<point> nodes;
    vector<scalar_t> weights;
    int deg;
    int sz;

    Quadrature volQuadratureRule;

    vector<scalar_t> refine_quad_x ;
    vector<scalar_t> refine_quad_y ;
    vector<scalar_t> refine_weight ;

    Matrix interpolate;
    Matrix nearMapping;
    Vector sqrtWeights;
    Vector legendreNorms;

    int numberOfSquares;
    int numberOfNodes;
    double dx;

};


#endif //ANISO2_GEOMETRY_H
