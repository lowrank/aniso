/*
 *
 * 1d and 2d (triangle) quadrature rule subroutine.
 *
 * copyright@ Yimin Zhong. yzhong@math.utexas.edu. All Rights Reserved.
 *
 */
#ifndef MODES_H
#define MODES_H


#include <vector>
#include <quadmath.h>
#include "bbfmm/utils.h"


using std::vector;

class Quadrature3 {
public:
    vector<scalar_t> points_x;
    vector<scalar_t> points_y;
    vector<scalar_t> weights;
    size_t degree;

    void resize(size_t n) {
        points_x.resize(n);
        points_y.resize(n);
        weights.resize(n);
    }
};


class Quadrature {
public:
    vector<scalar_t> points_x;
    vector<scalar_t> weights;

    void resize(size_t n) {
        points_x.resize(n);
        weights.resize(n);
    }
};

void get_legendre_data(size_t deg, Quadrature &table);

void affine(Quadrature &table);

void get_vr_data(size_t deg, Quadrature3 &table);

void affine(Quadrature3 &table);


#endif //ISORTES_MODES_H
