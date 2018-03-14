//
// Created by lurker on 2/4/18.
//

#include "Aniso.h"

Aniso::Aniso(int geometry_size, int geometry_degree,
             int kernel_size, double kernel_g,
             int kernel_degree, int _np, int _maxLevel):KernelFactory(geometry_size, geometry_degree,
                                              kernel_size, kernel_g, kernel_degree){
    np = _np;
    maxLevel = _maxLevel;
}

Aniso::~Aniso() {

}

void Aniso::displayKernelCacheSize() {
    for (int i = 0; i < kernelSize; ++i) {
        auto len = realParts[i].Cache.size();
        scalar_t cache = 0.; // avoid overflow
        for (int j = 0; j < len; ++j) {
            for (int k = 0 ; k < realParts[i].Cache[j].size(); ++k) {
                cache += realParts[i].Cache[j][k].row() * realParts[i].Cache[j][k].col() * sizeof(scalar_t);

            }
        }

        len = imagParts[i].Cache.size();
        for (int j = 0; j < len; ++j) {
            for (int k = 0 ; k < imagParts[i].Cache[j].size(); ++k) {
                cache += imagParts[i].Cache[j][k].row() * imagParts[i].Cache[j][k].col() * sizeof(scalar_t);

            }
        }


        auto res =  double(cache) / 1024.0 / 1024.0;
        if (res < 1024) {
            std::cout << res << " MByte" <<std::endl;
        }
        else {
            std::cout << res  / 1024.0 << " GByte" << std::endl;
        }
    }
}



