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



