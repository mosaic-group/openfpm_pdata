//
// Created by landfried on 06.12.21.
//

#include "Test1.hpp"
#include <iostream>

constexpr int dimension = 2;

typedef Test1<dimension> ParticleMethodType;
typedef InitialCondition1<ParticleMethodType> InitialConditionType;

int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    ParticleData<ParticleMethodType, InitialConditionType> particleData;

    TransitionCellList<ParticleMethodType, InitialConditionType> transition(particleData);

    while (!transition.stop(particleData)) {
        transition.run_step(particleData);
    }

    openfpm_finalize();
}