//
// Created by lurker on 2/4/18.
//

#include "Geometry.h"
/// Contructor of Geometry
///
/// @param size number of square per direction.
/// @param degree quadrature rule's degree
Geometry::Geometry(int size, int degree) {
    // get quadrature rule in 1D on [-1,1].
    get_legendre_data((size_t) degree, volQuadratureRule);

    deg = degree;
    dx = 1.0 / size;
    sz = size;
    numberOfSquares  = SQR(size);
    numberOfNodes = numberOfSquares * SQR(degree);

    nodes.resize((unsigned long) (numberOfSquares * SQR(degree)));
    weights.resize((unsigned long) (numberOfSquares * SQR(degree)));

    vector<scalar_t > quadrature_rule_x((unsigned long) SQR(degree));
    vector<scalar_t > quadrature_rule_y((unsigned long) SQR(degree));
    vector<scalar_t > quadrature_rule_w((unsigned long) SQR(degree));
    sqrtWeights.resize(SQR(degree));

    for (int r = 0; r < degree; ++r) {
        for (int c = 0; c < degree; ++c) {
            quadrature_rule_x[r * degree + c] = volQuadratureRule.points_x[r];
            quadrature_rule_y[r * degree + c] = volQuadratureRule.points_x[c];
            quadrature_rule_w[r * degree + c] = volQuadratureRule.weights[r] * volQuadratureRule.weights[c];
            sqrtWeights(r * degree + c) = sqrt(quadrature_rule_w[r * degree + c]);
        }
    }

    int nRefineLevel = 2;              // default value, can be adaptive.
    legendreNorms.resize(SQR(degree));

    /*
     * |----|----|----|
     * |    |    |    |  1
     * |----|----|----|     j
     * |    |    |    |  0
     * |----|----|----|
     *   0    1    2
     *        i
     */

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            for (int k = 0; k < SQR(degree); ++k){
                nodes[(i * size + j) * SQR(degree) + k].x =
                        (0.5 + i) * dx + 0.5 * (quadrature_rule_x[k]) * dx;
                nodes[(i * size + j) * SQR(degree) + k].y =
                        (0.5 + j) * dx + 0.5 * quadrature_rule_y[k] * dx;
                weights[(i * size + j) * SQR(degree) + k] =
                        quadrature_rule_w[k] * 0.25 * SQR(dx);
            }
        }
    }

    int coarseQuadratureSize = (int)SQR(degree);
    int refineQuadratureSize = coarseQuadratureSize;

    interpolate.resize(coarseQuadratureSize, coarseQuadratureSize);
    makeLegendreMatrix(interpolate, degree, quadrature_rule_x , quadrature_rule_y, quadrature_rule_w);

    vector<scalar_t> refine_quad_x = quadrature_rule_x;
    vector<scalar_t> refine_quad_y = quadrature_rule_y;
    vector<scalar_t> refine_weight = quadrature_rule_w;

    vector<scalar_t> refine_quad_x_tmp, refine_quad_y_tmp, refine_weight_tmp;

    for (int level = 0; level < nRefineLevel; ++level) {
        for (int id = 0; id < refineQuadratureSize; ++id) {
            refine_quad_x_tmp.push_back((refine_quad_x[id] + 1) / 2.0);
            refine_quad_y_tmp.push_back((refine_quad_y[id] + 1)/ 2.0);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);


            refine_quad_x_tmp.push_back((refine_quad_x[id] + 1) / 2.0);
            refine_quad_y_tmp.push_back((refine_quad_y[id] - 1)/ 2.0);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);

            refine_quad_x_tmp.push_back((refine_quad_x[id] - 1) / 2.0);
            refine_quad_y_tmp.push_back((refine_quad_y[id] + 1)/ 2.0);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);

            refine_quad_x_tmp.push_back((refine_quad_x[id] - 1) / 2.0);
            refine_quad_y_tmp.push_back((refine_quad_y[id] - 1)/ 2.0);
            refine_weight_tmp.push_back(refine_weight[id] / 4.0);
        }

        refine_quad_x = refine_quad_x_tmp;
        refine_quad_y = refine_quad_y_tmp;
        refine_weight = refine_weight_tmp;
        refineQuadratureSize *= 4;

        refine_quad_x_tmp.clear();
        refine_quad_y_tmp.clear();
        refine_weight_tmp.clear();
    }

    Matrix refinements(coarseQuadratureSize, refineQuadratureSize);      /* The interpolation on refined level */
    makeLegendreMatrix(refinements, degree, refine_quad_x, refine_quad_y, refine_weight);

    nearMapping.resize(refineQuadratureSize, coarseQuadratureSize);
    t_dgemm(1.0, refinements, interpolate, 0., nearMapping);
}

Geometry::~Geometry() {
#ifdef DEBUG
    std::cout << "Geometry detached." << std::endl;
#endif
}

/// makeLegendreMatrix: construct the Legendre matrix with quadrature points and weights.
///
/// @param K output matrix
/// @param N number of degree
/// @param x x coordinate
/// @param y y coordiante
/// @param w weights
void Geometry::makeLegendreMatrix(Matrix  &K, int N, vector<scalar_t> &x, vector<scalar_t> &y, vector<scalar_t> &w){
    int row = 0;
    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < N; ++k) {
            for (int I = 0; I < x.size(); ++I) {
                K(row, I) = (scalar_t)(legendre(n, x[I]) * legendre(k, y[I]) * sqrt(w[I]));
            }
            ++row;
        }
    }

    for (row = 0; row < K.row(); ++row) {
        scalar_t nrm = 0.;
        for (int col = 0; col < K.col(); ++col) {
            nrm += SQR(K(row, col));
        }
        nrm = sqrt(nrm);

        legendreNorms(row) = nrm;

        assert(nrm > EPS);
        for (int col = 0; col < K.col(); ++col) {
            K(row, col) /= nrm;
        }
    }
}


