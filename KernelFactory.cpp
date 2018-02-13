//
// Created by lurker on 2/4/18.
//

#include "KernelFactory.h"

KernelFactory::KernelFactory(int geometry_size, int geometry_degree, int kernel_size,
                             double kernel_g, int kernel_degree):Geometry(geometry_size, geometry_degree) {
    assert(size >= 1);
    realParts.resize((unsigned long) (2 * kernel_size - 1));
    kernelSize = 2 * kernel_size - 1;

    get_legendre_data((size_t) kernel_degree, singQuadratureRule);
    get_legendre_data((size_t) geometry_degree, lineQuadratureRule);


    // set sigma_s, sigma_t.
    sigma_s.resize(numberOfNodes);
    sigma_t.resize(numberOfNodes);

    sigma_s_coeff.resize((unsigned long) numberOfSquares);
    sigma_t_coeff.resize((unsigned long) numberOfSquares);

    for (int i = 0; i < numberOfSquares; ++i) {
        sigma_s_coeff[i].resize(SQR(deg));
        sigma_t_coeff[i].resize(SQR(deg));
    }

    np = 4;
    maxLevel = 10;
    maxPoint = 16;

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
        auto row0 = getCol(y0);
        double ybar = double(row0)/side;
        double xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
        return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
    }
    else if ((row0 == row1 - 1) && (col0 == col1)) {
        auto row1 = getCol(y1);
        double ybar = double(row1)/side;
        double xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
        return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 + 1) && (row0 == row1)) {
        auto col0 = getCol(x0);
        double xbar = double(col0)/side;
        double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
        return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
    }
    else if ((col0 == col1 - 1) && (row0 == row1)) {
        auto col1 = getCol(x1);
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
                load(n * deg + k) = legendre(n, x) * legendre(k, y) / legendreNorms(n * deg + k);
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

/// helper function to calculate local line integral.
/// \param p
/// \param q
/// \return
scalar_t KernelFactory::integral_helper(point &p, point &q) {
    return integral_helper(p.x, p.y, q.x, q.y);
}

void KernelFactory::makeKernels() {
    for (int i = 0; i < kernelSize; ++i) {
        // capture this and the kernel index
        realParts[i].eval = [this,i](point& a, point& b) {
            scalar_t dist = sqrt(SQR(a.x - b.x) + SQR(a.y - b.y));
            scalar_t ang  = atan2(a.y - b.y, a.x - b.x);
            return dist == 0. ? 0. : exp(-lineIntegral(a, b))*cos(i * ang)/dist;
        };
    }

}

void KernelFactory::runKernels(Vector& f) {
    // not finished yet.
    //todo: f is values, assign weights.
    for (int i = 0; i < kernelSize; ++i) {
        realParts[i].initialize(np, nodes, nodes, f,
                                (index_t) nodes.size(),
                                (index_t) nodes.size(), np * np, maxLevel);
        Vector ret;
        realParts[i].run(ret);
    }
}

int  KernelFactory::getRow(double y)  {
    return int(floor(y * sz));
}

int  KernelFactory::getCol(double x)  {
    return int(floor(x * sz));
}