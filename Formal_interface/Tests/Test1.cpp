//
// Created by landfried on 06.12.21.
//

#include "Test1.hpp"
#include <iostream>
#include "../InitialCondition.hpp"
constexpr int dimension = 2;

typedef Test1<dimension> ParticleMethodType;
typedef SimulationParams1<ParticleMethodType> SimulationParametersType;

int main(int argc, char* argv[]) {

//    typedef InitialConditionMesh ictype;
//    InitialCondition_Impl<ictype> IC_Impl;
//    IC_Impl.initialization();


    openfpm_init(&argc,&argv);

    ParticleData<ParticleMethodType, SimulationParametersType> particleData;

    TransitionCellList<ParticleMethodType, SimulationParametersType> transition(particleData);

    while (!transition.stop(particleData)) {
        transition.run_step(particleData);
    }

    openfpm_finalize();
}