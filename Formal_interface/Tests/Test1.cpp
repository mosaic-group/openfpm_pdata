//
// Created by landfried on 06.12.21.
//

#include "Test1.hpp"
#include <iostream>

constexpr int dimension = 2;

typedef Test1<dimension> ParticleMethodType;

int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

//    PD<ParticleMethodType> pd;

    ParticleData<ParticleMethodType/*, InitializationType*/> particleData;

    TransitionCellList<ParticleMethodType> transition(particleData);

    while (!transition.stop(particleData)) {
        transition.run_step(particleData);
    }

    openfpm_finalize();
}