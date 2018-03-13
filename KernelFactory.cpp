//
// Created by lurker on 2/4/18.
//

#include "KernelFactory.h"

KernelFactory::KernelFactory(int geometry_size, int geometry_degree, int kernel_size,
                             double kernel_g, int kernel_degree):Geometry(geometry_size, geometry_degree) {
    assert(size >= 1);
    realParts.resize((unsigned long) (2 * kernel_size - 1));
    anisotropy.resize(kernel_size);
    kernelSize = 2 * kernel_size - 1;

    get_legendre_data((size_t) kernel_degree, singQuadratureRule); // Duffy
    affine(singQuadratureRule);

    for (int i = 0; i < kernel_size; ++i) {
        anisotropy(i) = (pow(kernel_g, i) - pow(kernel_g, kernel_size))/(1 - pow(kernel_g, kernel_size));
    }

    // set sigma_s, sigma_t.
    sigma_s.resize(numberOfNodes);
    sigma_t.resize(numberOfNodes);

    sigma_s_Coeff.resize((unsigned long) numberOfSquares);
    sigma_t_Coeff.resize((unsigned long) numberOfSquares);

    for (int i = 0; i < numberOfSquares; ++i) {
        sigma_s_Coeff[i].resize(SQR(deg));
        sigma_t_Coeff[i].resize(SQR(deg));
    }



    nearInteractions.resize((unsigned long) kernelSize);
    singInteractions.resize((unsigned long) kernelSize);
    for (int kernelId = 0; kernelId < kernelSize; ++kernelId) {
        nearInteractions[kernelId].resize((unsigned long) (numberOfSquares * SQR(deg)));
        singInteractions[kernelId].resize((unsigned long) (numberOfSquares * SQR(deg)));
        for (int nodeId = 0; nodeId < numberOfSquares * SQR(deg); ++nodeId) {
            nearInteractions[kernelId][nodeId].resize(3);
            singInteractions[kernelId][nodeId].resize(8 * SQR(singQuadratureRule.weights.size()));
            for (int rowId = 0; rowId < 3; ++rowId) {
                nearInteractions[kernelId][nodeId][rowId].resize(3);
                for (int colId = 0; colId < 3; ++colId) {
                    nearInteractions[kernelId][nodeId][rowId][colId].resize((unsigned long) refineQuadratureSize);
                }
            }
        }
    }


}


KernelFactory::~KernelFactory() {

}

/// lineIntegral
/// \param x0
/// \param y0
/// \param x1
/// \param y1
/// \return
scalar_t KernelFactory::lineIntegral(double x0, double y0, double x1, double y1) {
    auto col0 = getRow(x0);
    auto col1 = getRow(x1);
    auto row0 = getCol(y0);
    auto row1 = getCol(y1);
    auto side = 1.0/dx;
    // 9 cases
    if ((row0 == row1) && (col0 == col1)) {
        return integral_helper(x0, y0, x1, y1);
    }
    else if ((row0 == row1 + 1) && (col0 == col1)) {
        double ybar = double(row0)/side;
        double xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
        return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
    }
    else if ((row0 == row1 - 1) && (col0 == col1)) {
        double ybar = double(row1)/side;
        double xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
        return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 + 1) && (row0 == row1)) {
        double xbar = double(col0)/side;
        double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 - 1) && (row0 == row1)) {
        double xbar = double(col1)/side;
        double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 + 1) && (row0 == row1 + 1)) {
        double xbar = double(col0)/side;
        double ybar2 = double(row0)/side;
        double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar < xbar2) {
            return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
        }
        else {
            return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);

        }
    }
    else if ((col0 == col1 + 1) && (row0 == row1 - 1)) {
        double xbar = double(col0)/side;
        double ybar2 = double(row1)/side;
        double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar < xbar2) {
            return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
        }
        else {
            return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);

        }

    }
    else if ((col0 == col1 - 1) && (row0 == row1 + 1)) {
        double xbar = double(col1)/side;
        double ybar2 = double(row0)/side;
        double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar > xbar2) {
            return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
        }
        else {
            return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);
        }
    }
    else if ((col0 == col1 - 1) && (row0 == row1 - 1)) {
        double xbar = double(col1)/side;
        double ybar2 = double(row1)/side;
        double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

        if (xbar > xbar2) {
            return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
        }
        else {
            return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);

        }
    }
    else {
        double xm = (x0 + x1)/2;
        double ym = (y0 + y1)/2;
        return lineIntegral(x0, y0, xm, ym) + lineIntegral(xm, ym, x1, y1);
    }
}

/// integral_helper
/// \param x0
/// \param y0
/// \param x1
/// \param y1
/// \return
scalar_t KernelFactory::integral_helper(double x0, double y0, double x1, double y1) {
    Vector load(SQR(deg));
    int col = getRow((x0+x1)/2);
    int row = getCol((y0+y1)/2);
    scalar_t ret = 0.;
    for(int i =  0; i < deg; ++i) {
        scalar_t x = (x0 + x1) / 2 + (x0 - x1) / 2 * volQuadratureRule.points_x[i];
        scalar_t y = (y0 + y1) / 2 + (y0 - y1) / 2 * volQuadratureRule.points_x[i];
        for (int n = 0; n < deg; ++n) {
            for (int k = 0; k < deg; ++k) {
                load(n * deg + k) = legendre((unsigned int) n, x) * legendre((unsigned int) k, y) / legendreNorms(n * deg + k);
            }
        }
        ret += ddot(load, sigma_t_Coeff[col * sz + row]) * volQuadratureRule.weights[i];
    }
    return ret * sqrt(SQR(x0 - x1) + SQR(y0 - y1))/2.0;
}

void KernelFactory::interpolation() {

#ifdef RUN_OMP
#pragma omp parallel for schedule(static, CHUNKSIZE)  num_threads(omp_get_max_threads())
#endif
    for (int i = 0; i < numberOfSquares;++i) {
        Vector load_t(SQR(deg));
        Vector load_s(SQR(deg));
        for (int j = 0; j < SQR(deg);++j) {
            load_t(j) = sqrtWeights(j) * sigma_t(i * SQR(deg) + j);
            load_s(j) = sqrtWeights(j) * sigma_s(i * SQR(deg) + j);
        }
        dgemv(1.0, interpolate, load_t, 0, sigma_t_Coeff[i]);
        dgemv(1.0, interpolate, load_s, 0, sigma_s_Coeff[i]);
    }
}

/// must cache such scalar values, since it is very large cost for many particle interactions.
/// \param p
/// \param q
/// \return
scalar_t KernelFactory::lineIntegral(point &p, point &q) {
    return lineIntegral(p.x, p.y, q.x, q.y);
}

/// generate the kernels.
/// \param
/// \return
void KernelFactory::makeKernels() {
    for (int i = 0; i < kernelSize; ++i) {
        // capture this and the kernel index
        realParts[i].eval = [this,i](point& a, point& b) {
            scalar_t dist = sqrt(SQR(a.x - b.x) + SQR(a.y - b.y));
            scalar_t ang  = atan2(a.y - b.y, a.x - b.x);
            // compute the traditional interaction by fmm.
            return dist == 0. ? 0. :  exp(-lineIntegral(a, b))*cos(i * ang)/dist;
        };
    }
}


/// \param f
void KernelFactory::runKernels(Vector& f, Vector& ret) {
    // not finished yet.
    for (int i = 0; i < kernelSize; ++i) {
        realParts[i].initialize(np, nodes, nodes, f,
                                (index_t) nodes.size(),
                                (index_t) nodes.size(), np * np, maxLevel);
        realParts[i].run(ret);
    }
}
/// input can be 0 vector only to cache the kernels.
/// f can be dropped.
/// \param f
void KernelFactory::runKernelsCache(Vector& f, Vector& ret) {
    // not finished yet.
    for (int i = 0; i < kernelSize; ++i) {
        realParts[i].initialize(np, nodes, nodes, f,
                                (index_t) nodes.size(),
                                (index_t) nodes.size(), np * np, maxLevel);

        realParts[i].runCache(ret);
    }
}

/// Run the kernels with rhs sources. Should be something like Toeplitz matrix product.
/// \param f
void KernelFactory::runKernelsFast(Vector& f, Vector& ret) {
    // not finished yet.
    for (int i = 0; i < kernelSize; ++i) {
        realParts[i].initialize(np, nodes, nodes, f,
                                (index_t) nodes.size(),
                                (index_t) nodes.size(), np * np, maxLevel);


        realParts[i].runFast(ret);
    }
}
/// get the row number for y coordinate
/// \param y
/// \return
int  KernelFactory::getRow(double y)  {
    return int(floor(y * sz));
}

/// get the col number for x coordinate
/// \param x
/// \return
int  KernelFactory::getCol(double x)  {
    return int(floor(x * sz));
}

/// costy, remove all nearby interactions.
/// \param f
void KernelFactory::nearRemoval(Vector &f, Vector& ret) {
    for (int i = 0; i < kernelSize; ++i) {
#ifdef RUN_OMP
#pragma omp parallel for schedule(static, CHUNKSIZE) collapse(2) num_threads(omp_get_max_threads())
#endif
        for (int targetSquareId = 0; targetSquareId < numberOfSquares; ++targetSquareId) {
            for (int targetQuadratureId = 0; targetQuadratureId < SQR(deg); ++targetQuadratureId) {

                int targetId = targetSquareId * (SQR(deg)) + targetQuadratureId;
                int targetSquareRow = targetSquareId / sz;
                int targetSquareCol = targetSquareId - sz * targetSquareRow;

                for (int nearSourceSquareRowId = -1; nearSourceSquareRowId < 2; ++nearSourceSquareRowId) {
                    for (int nearSourceSquareColId = -1; nearSourceSquareColId < 2; ++nearSourceSquareColId) {
                        int nearSourceSquareId = targetSquareId + nearSourceSquareRowId * sz + nearSourceSquareColId;

//                        if (nearSourceSquareId == targetSquareId) continue;

                        if (nearSourceSquareRowId + targetSquareRow >= 0 && nearSourceSquareRowId + targetSquareRow < sz) {
                            if (nearSourceSquareColId + targetSquareCol >= 0 && nearSourceSquareColId + targetSquareCol < sz) {
                                for (int nearSourceQuadratureId = 0;
                                     nearSourceQuadratureId < SQR(deg); ++nearSourceQuadratureId) {
                                    int nearSourceId = nearSourceSquareId * (SQR(deg)) + nearSourceQuadratureId;
                                    ret(targetId) -=
                                            realParts[i].eval(nodes[nearSourceId], nodes[targetId]) * f(nearSourceId);

                                }
                            }
                        }
                    }
                }

            }
        }
    }
}

/// costy, add refined interactions.
/// \param f
void KernelFactory::refineAddOnCache(Vector &f, Vector& ret) {
    for (int i = 0; i < kernelSize; ++i) {
#ifdef RUN_OMP
#pragma omp parallel for schedule(static, CHUNKSIZE) collapse(2) num_threads(omp_get_max_threads())
#endif
        for (int targetSquareId = 0; targetSquareId < numberOfSquares; ++targetSquareId) {
            for (int targetQuadratureId = 0; targetQuadratureId < SQR(deg); ++targetQuadratureId) {

                int targetId = targetSquareId * (SQR(deg)) + targetQuadratureId;
                int targetSquareRow = targetSquareId / sz;
                int targetSquareCol = targetSquareId - sz * targetSquareRow;

                for (int nearSourceSquareRowId = -1; nearSourceSquareRowId < 2; ++nearSourceSquareRowId) {
                    for (int nearSourceSquareColId = -1; nearSourceSquareColId < 2; ++nearSourceSquareColId) {
                        int nearSourceSquareId = targetSquareId + nearSourceSquareRowId * sz + nearSourceSquareColId;

                        if (nearSourceSquareId == targetSquareId) continue;

                        if (nearSourceSquareRowId + targetSquareRow >= 0 &&
                            nearSourceSquareRowId + targetSquareRow < sz) {
                            if (nearSourceSquareColId + targetSquareCol >= 0 &&
                                nearSourceSquareColId + targetSquareCol < sz) {

                                Vector oldValues(coarseQuadratureSize);
                                Vector newValues(refineQuadratureSize);

                                for (int sourceQuadratureId = 0; sourceQuadratureId < coarseQuadratureSize; ++sourceQuadratureId) {
                                    oldValues(sourceQuadratureId) =
                                            f(nearSourceSquareId * coarseQuadratureSize +
                                                      sourceQuadratureId) / sqrtWeights(sourceQuadratureId);
                                }

                                dgemv(1.0, nearMapping, oldValues, 0., newValues);

                                for (int nearSourceQuadratureId = 0;
                                     nearSourceQuadratureId < refineQuadratureSize; ++nearSourceQuadratureId){

                                    scalar_t lambda = (scalar_t) refine_quad_x[nearSourceQuadratureId];
                                    scalar_t mu     = (scalar_t) refine_quad_y[nearSourceQuadratureId];
                                    scalar_t w      = (scalar_t) refine_weight[nearSourceQuadratureId];

                                    point cur_nearSourcePoint = {
                                            (0.5 + (targetSquareRow + nearSourceSquareRowId)) * dx + 0.5 * (lambda) * dx,
                                            (0.5 + (targetSquareCol + nearSourceSquareColId )) * dx + 0.5 * (mu) * dx
                                    };

                                    // Cache the nearRefinement interaction.
                                    nearInteractions[i][targetId][nearSourceSquareRowId+1][nearSourceSquareColId+1][nearSourceQuadratureId] =
                                            realParts[i].eval(cur_nearSourcePoint, nodes[targetId]) * sqrt(w);
                                    ret(targetId) += nearInteractions[i][targetId][nearSourceSquareRowId+1][nearSourceSquareColId+1][nearSourceQuadratureId]
                                                    * newValues(nearSourceQuadratureId);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


void KernelFactory::refineAddOnFast(Vector &f, Vector& ret) {
    for (int i = 0; i < kernelSize; ++i) {
#ifdef RUN_OMP
#pragma omp parallel for schedule(static, CHUNKSIZE) collapse(2) num_threads(omp_get_max_threads())
#endif
        for (int targetSquareId = 0; targetSquareId < numberOfSquares; ++targetSquareId) {
            for (int targetQuadratureId = 0; targetQuadratureId < SQR(deg); ++targetQuadratureId) {

                int targetId = targetSquareId * (SQR(deg)) + targetQuadratureId;
                int targetSquareRow = targetSquareId / sz;
                int targetSquareCol = targetSquareId - sz * targetSquareRow;

                for (int nearSourceSquareRowId = -1; nearSourceSquareRowId < 2; ++nearSourceSquareRowId) {
                    for (int nearSourceSquareColId = -1; nearSourceSquareColId < 2; ++nearSourceSquareColId) {
                        int nearSourceSquareId = targetSquareId + nearSourceSquareRowId * sz + nearSourceSquareColId;

                        if (nearSourceSquareId == targetSquareId) continue;

                        if (nearSourceSquareRowId + targetSquareRow >= 0 &&
                            nearSourceSquareRowId + targetSquareRow < sz) {
                            if (nearSourceSquareColId + targetSquareCol >= 0 &&
                                nearSourceSquareColId + targetSquareCol < sz) {

                                Vector oldValues(coarseQuadratureSize);
                                Vector newValues(refineQuadratureSize);

                                for (int sourceQuadratureId = 0; sourceQuadratureId < coarseQuadratureSize; ++sourceQuadratureId) {
                                    oldValues(sourceQuadratureId) =
                                            f(nearSourceSquareId * coarseQuadratureSize +
                                              sourceQuadratureId) / sqrtWeights(sourceQuadratureId);
                                }

                                dgemv(1.0, nearMapping, oldValues, 0., newValues);

                                for (int nearSourceQuadratureId = 0;
                                     nearSourceQuadratureId < refineQuadratureSize; ++nearSourceQuadratureId){

                                    ret(targetId) += nearInteractions[i][targetId][nearSourceSquareRowId+1][nearSourceSquareColId+1][nearSourceQuadratureId]
                                                     * newValues(nearSourceQuadratureId);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void KernelFactory::singularAddCache(vector<Vector>& f_coeff, Vector& ret) {
    for (int i = 0; i < kernelSize; ++i) {
        // no precomputation is acceptable.
#ifdef RUN_OMP
#pragma omp parallel for schedule(static, CHUNKSIZE) collapse(2)  num_threads(omp_get_max_threads())
#endif
        for (int targetSquareId = 0; targetSquareId < numberOfSquares; ++targetSquareId) {
            for (int targetQuadratureId = 0; targetQuadratureId < SQR(deg); ++targetQuadratureId) {

                int col = targetSquareId / sz;
                int row = targetSquareId - col * sz;

                int targetId = targetSquareId * (SQR(deg)) + targetQuadratureId;

                for (int sourceQuadratureId = 0; sourceQuadratureId < 8 * SQR(singQuadratureRule.weights.size()); ++sourceQuadratureId) {
                    // scaled coordinates
                    scalar_t x =  (0.5 + col) * dx + 0.5 * (singX[targetQuadratureId][sourceQuadratureId]) * dx;
                    scalar_t y =  (0.5 + row) * dx + 0.5 * (singY[targetQuadratureId][sourceQuadratureId]) * dx;
                    scalar_t w = singW[targetQuadratureId][sourceQuadratureId] * SQR(dx) / 4.0;

                    point cur_point = {x, y};

                    Vector load(SQR(deg));
                    for (int n = 0; n < deg; ++n) {
                        for (int k = 0; k < deg; ++k) {
                            load(n * deg + k) = legendre((unsigned int) n, x) * legendre((unsigned int) k, y) / legendreNorms(n * deg + k);
                        }
                    }

                    singInteractions[i][targetId][sourceQuadratureId] = realParts[i].eval(cur_point, nodes[targetId]) * w;
                    ret(targetId) += ddot(load, f_coeff[col * sz + row]) *
                            singInteractions[i][targetId][sourceQuadratureId];

                }
            }
        }
    }
}

void KernelFactory::singularAddFast(vector<Vector>& f_coeff, Vector& ret) {
    for (int i = 0; i < kernelSize; ++i) {
        // no precomputation is acceptable.
#ifdef RUN_OMP
#pragma omp parallel for schedule(static, CHUNKSIZE) collapse(2)  num_threads(omp_get_max_threads())
#endif
        for (int targetSquareId = 0; targetSquareId < numberOfSquares; ++targetSquareId) {
            for (int targetQuadratureId = 0; targetQuadratureId < SQR(deg); ++targetQuadratureId) {

                int col = targetSquareId / sz;
                int row = targetSquareId - col * sz;
                int targetId = targetSquareId * (SQR(deg)) + targetQuadratureId;

                for (int sourceQuadratureId = 0; sourceQuadratureId < 8 * SQR(singQuadratureRule.weights.size()); ++sourceQuadratureId) {
                    // scaled coordinates
                    scalar_t x =  (0.5 + col) * dx + 0.5 * (singX[targetQuadratureId][sourceQuadratureId]) * dx;
                    scalar_t y =  (0.5 + row) * dx + 0.5 * (singY[targetQuadratureId][sourceQuadratureId]) * dx;

                    point cur_point = {x, y};

                    Vector load(SQR(deg));
                    for (int n = 0; n < deg; ++n) {
                        for (int k = 0; k < deg; ++k) {
                            load(n * deg + k) = legendre((unsigned int) n, x) * legendre((unsigned int) k, y) / legendreNorms(n * deg + k);
                        }
                    }

                    ret(targetId) += ddot(load, f_coeff[col * sz + row]) *
                                     singInteractions[i][targetId][sourceQuadratureId];

                }
            }
        }

        //std::cout << std::setprecision(16) << ret(0) << " " << ret((int) (1)) << std::endl;
    }
}


void KernelFactory::duffy_transform(vector<scalar_t> points, vector<scalar_t> &X, vector<scalar_t> &Y, vector<scalar_t> &W) {

    scalar_t a = points[0];
    scalar_t b = points[1];

    scalar_t a11 = points[2] - points[0];
    scalar_t a12 = points[4] - points[2];
    scalar_t a21 = points[3] - points[1];
    scalar_t a22 = points[5] - points[3];

    Matrix A(2, 2);
    A(0, 0) = a11;
    A(0, 1) = a12;
    A(1, 0) = a21;
    A(1, 1) = a22;

    scalar_t detA = a11 * a22 - a12 * a21;

    size_t nsp =singQuadratureRule.points_x.size();

    vector<scalar_t> sx(nsp * nsp);
    vector<scalar_t> sy(nsp * nsp);
    vector<scalar_t> sw(nsp * nsp);


    int id = 0;
    for (int ri = 0; ri < nsp; ++ri) {
        for (int ci = 0; ci < nsp; ++ci) {
            sx[id] = singQuadratureRule.points_x[ri];
            sy[id] = singQuadratureRule.points_x[ci];
            sw[id++] = singQuadratureRule.weights[ri] * singQuadratureRule.weights[ci];
        }
    }

    vector<scalar_t> Z1x(nsp * nsp);
    vector<scalar_t> Z1y(nsp * nsp);
    vector<scalar_t> Z1w(nsp * nsp);

    X.resize(nsp * nsp);
    Y.resize(nsp * nsp);
    W.resize(nsp * nsp);


    for (int i = 0; i < nsp * nsp; ++i) {
        scalar_t u = sx[i];
        scalar_t v = sy[i];
        scalar_t w = sw[i];

        Z1x[i] = u;
        Z1y[i] = u * v;
        Z1w[i] = w * u;

        u = Z1x[i];
        v = Z1y[i];
        w = Z1w[i];

        scalar_t x = a11 * u + a12 * v + a;
        scalar_t y = a21 * u + a22 * v + b;
        scalar_t wr = detA * w;

        X[i] = x;
        Y[i] = y;
        W[i] = wr;
    }
}

void KernelFactory::singPrecompute() {
    vector<vector< vector<scalar_t > > > X, Y, W;

    X.resize((unsigned long) SQR(deg));
    Y.resize((unsigned long) SQR(deg));
    W.resize((unsigned long) SQR(deg));

    for (int targetId = 0; targetId < SQR(deg); ++targetId) {
        X[targetId].resize(8); Y[targetId].resize(8); W[targetId].resize(8);
    }

    for (int targetQuadratureRowId = 0; targetQuadratureRowId < deg; ++targetQuadratureRowId) {
        for (int targetQuadratureColId = 0; targetQuadratureColId < deg; ++targetQuadratureColId) {

            int targetId = targetQuadratureRowId * deg + targetQuadratureColId;

            scalar_t x = volQuadratureRule.points_x[targetQuadratureRowId];
            scalar_t y = volQuadratureRule.points_x[targetQuadratureColId];

            // 8 triangles.
            // 1st. x,y -- 1, y -- 1,1
            // 2nd, x,y -- 1,1 -- x, 1
            // 3rd, x,y -- x, 1 -- -1,1
            // 4th, x,y -- -1,1 -- -1, y
            // 5th, x,y -- -1,y -- -1,-1
            // 6th, x,y -- -1,-1 -- x, -1
            // 7th, x,y -- x, -1 -- 1,-1
            // 8th, x,y -- 1,-1 -- 1, y

            duffy_transform({x, y, 1. , y  , 1. ,   1.}, X[targetId][0], Y[targetId][0], W[targetId][0]);
            duffy_transform({x, y, 1. , 1  , x  ,   1.}, X[targetId][1], Y[targetId][1], W[targetId][1]);
            duffy_transform({x, y, x  , 1. , -1.,   1.}, X[targetId][2], Y[targetId][2], W[targetId][2]);
            duffy_transform({x, y, -1., 1. , -1.,    y}, X[targetId][3], Y[targetId][3], W[targetId][3]);
            duffy_transform({x, y, -1., y  , -1.,  -1.}, X[targetId][4], Y[targetId][4], W[targetId][4]);
            duffy_transform({x, y, -1., -1 , x  ,  -1.}, X[targetId][5], Y[targetId][5], W[targetId][5]);
            duffy_transform({x, y, x  , -1., 1. ,  -1.}, X[targetId][6], Y[targetId][6], W[targetId][6]);
            duffy_transform({x, y, 1. , -1 , 1. ,    y}, X[targetId][7], Y[targetId][7], W[targetId][7]);
        }
    }

    singX.resize((unsigned long) SQR(deg));
    singY.resize((unsigned long) SQR(deg));
    singW.resize((unsigned long) SQR(deg));

    for (int targetId = 0; targetId < SQR(deg); targetId++) {
        singX[targetId].resize(8 * SQR(singQuadratureRule.weights.size()));
        singY[targetId].resize(8 * SQR(singQuadratureRule.weights.size()));
        singW[targetId].resize(8 * SQR(singQuadratureRule.weights.size()));
        for (int triId = 0; triId < 8; ++triId) {
            for (int quadId = 0; quadId < SQR(singQuadratureRule.weights.size()); ++quadId) {
                singX[targetId][triId * SQR(singQuadratureRule.weights.size()) + quadId] = X[targetId][triId][quadId];
                singY[targetId][triId * SQR(singQuadratureRule.weights.size()) + quadId] = Y[targetId][triId][quadId];
                singW[targetId][triId * SQR(singQuadratureRule.weights.size()) + quadId] = W[targetId][triId][quadId];
            }
        }
    }

}

void KernelFactory::interpolation(Vector &h, vector<Vector> &h_coeff) {
    // coefficients
    h_coeff.resize((unsigned long) numberOfSquares);
    for (int i = 0; i < numberOfSquares; ++i) {
        h_coeff[i].resize(SQR(deg));
    }

#ifdef RUN_OMP
#pragma omp parallel for schedule(static, CHUNKSIZE)  num_threads(omp_get_max_threads())
#endif
    for (int i = 0; i < numberOfSquares;++i) {
        Vector h_t(SQR(deg));
        for (int j = 0; j < SQR(deg);++j) {
            h_t(j) = sqrtWeights(j) * h(i * SQR(deg) + j);
        }
        dgemv(1.0, interpolate, h_t, 0, h_coeff[i]);
    }
}

