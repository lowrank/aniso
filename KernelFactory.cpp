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
    get_legendre_data((size_t) geometry_degree, lineQuadratureRule); // Usual case

    for (int i = 0; i < kernel_size; ++i) {
        anisotropy(i) = (pow(kernel_g, i) - pow(kernel_g, kernel_size))/(1 - pow(kernel_g, kernel_size));
    }


    // set sigma_s, sigma_t.
    sigma_s.resize(numberOfNodes);
    sigma_t.resize(numberOfNodes);

    sigma_s_coeff.resize((unsigned long) numberOfSquares);
    sigma_t_coeff.resize((unsigned long) numberOfSquares);

    for (int i = 0; i < numberOfSquares; ++i) {
        sigma_s_coeff[i].resize(SQR(deg));
        sigma_t_coeff[i].resize(SQR(deg));
    }

    nearInteractions.resize((unsigned long) kernelSize);
    for (int kernelId = 0; kernelId < kernelSize; ++kernelId) {
        nearInteractions[kernelId].resize((unsigned long) (numberOfSquares * SQR(deg)));
        for (int nodeId = 0; nodeId < numberOfSquares * SQR(deg); ++nodeId) {
            nearInteractions[kernelId][nodeId].resize(3);
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
        ret += ddot(load, sigma_t_coeff[col * sz + row]) * volQuadratureRule.weights[i];
    }
    return ret * sqrt(SQR(x0 - x1) + SQR(y0 - y1))/2.0;
}

void KernelFactory::interpolation() {
    Vector load_t(SQR(deg));
    Vector load_s(SQR(deg));
#pragma omp parallel for
    for (int i = 0; i < numberOfSquares;++i) {
        for (int j = 0; j < SQR(deg);++j) {
            load_t(j) = sqrtWeights(j) * sigma_t(i * SQR(deg) + j);
            load_s(j) = sqrtWeights(j) * sigma_s(i * SQR(deg) + j);
        }
        dgemv(1.0, interpolate, load_t, 0, sigma_t_coeff[i]);
        dgemv(1.0, interpolate, load_s, 0, sigma_s_coeff[i]);
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
void KernelFactory::runKernels(Vector& f) {
    // not finished yet.
    for (int i = 0; i < kernelSize; ++i) {
        realParts[i].initialize(np, nodes, nodes, f,
                                (index_t) nodes.size(),
                                (index_t) nodes.size(), np * np, maxLevel);
        Vector ret;
        realParts[i].run(ret);
        std::cout << ret(0) << " " << ret((int) (nodes.size() / 2)) << std::endl;
    }
}
/// input can be 0 vector only to cache the kernels.
/// f can be dropped.
/// \param f
void KernelFactory::runKernelsCache(Vector& f) {
    // not finished yet.
    for (int i = 0; i < kernelSize; ++i) {
        realParts[i].initialize(np, nodes, nodes, f,
                                (index_t) nodes.size(),
                                (index_t) nodes.size(), np * np, maxLevel);
        Vector ret;
        realParts[i].runCache(ret);

        std::cout << ret(0) << " " << ret((int) (nodes.size() / 2)) << std::endl;
    }
}

/// Run the kernels with rhs sources. Should be something like Toeplitz matrix product.
/// \param f
void KernelFactory::runKernelsFast(Vector& f) {
    // not finished yet.
    for (int i = 0; i < kernelSize; ++i) {
        realParts[i].initialize(np, nodes, nodes, f,
                                (index_t) nodes.size(),
                                (index_t) nodes.size(), np * np, maxLevel);

        Vector ret;
        realParts[i].runFast(ret);
        std::cout << ret(0) << " " << ret((int) (nodes.size() / 2)) << std::endl;
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
void KernelFactory::nearRemoval(Vector &f) {
    for (int i = 0; i < kernelSize; ++i) {
        Vector ret((int) nodes.size());
        setValue(ret, 0.);
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
                                    ret(targetId) +=
                                            realParts[i].eval(nodes[nearSourceId], nodes[targetId]) * f(nearSourceId);

                                }
                            }
                        }
                    }
                }

            }
        }
        std::cout << ret(0) << " " << ret((int) (nodes.size() / 2)) << std::endl;
    }
}

/// costy, add refined interactions.
/// \param f
void KernelFactory::refineAddOnCache(Vector &f) {
    for (int i = 0; i < kernelSize; ++i) {
        Vector ret((int) nodes.size());
        setValue(ret, 0.);

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
                                    nearInteractions[i][targetId][nearSourceSquareRowId+1][nearSourceSquareColId+1][nearSourceQuadratureId] = realParts[i].eval(cur_nearSourcePoint, nodes[targetId]);
                                    ret(targetId) += nearInteractions[i][targetId][nearSourceSquareRowId+1][nearSourceSquareColId+1][nearSourceQuadratureId]
                                                     * sqrt(w) * newValues(nearSourceQuadratureId);

                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << std::setprecision(16) << ret(0) << " " << ret((int) (nodes.size() / 2)) << std::endl;
    }
}


void KernelFactory::refineAddOnFast(Vector &f) {
    for (int i = 0; i < kernelSize; ++i) {
        Vector ret((int) nodes.size());
        setValue(ret, 0.);

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
                                    //nearInteractions[i][targetId][nearSourceSquareRowId][nearSourceSquareColId][nearSourceQuadratureId] = realParts[i].eval(cur_nearSourcePoint, nodes[targetId]);
                                    ret(targetId) += nearInteractions[i][targetId][nearSourceSquareRowId+1][nearSourceSquareColId+1][nearSourceQuadratureId]
                                                     * sqrt(w) * newValues(nearSourceQuadratureId);

                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << std::setprecision(16) << ret(0) << " " << ret((int) (nodes.size() / 2)) << std::endl;
    }
}

void KernelFactory::singularAdd(Vector &f) {

}
