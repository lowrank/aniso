//
// Created by lurker on 2/4/18.
//

#ifndef ANISO2_ANISO_H
#define ANISO2_ANISO_H

#include "Geometry.h"
#include "KernelFactory.h"
#include "utility/Profiler.h"

class Aniso : public KernelFactory {
public:
    Aniso(int geometry_size, int geometry_degree,
          int kernel_size, double kernel_g, int kernel_degree, int np, int maxLevel);
    ~Aniso();


    void CacheSize();



};


#endif //ANISO2_ANISO_H
