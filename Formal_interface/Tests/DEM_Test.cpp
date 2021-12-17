//
// Created by landfried on 06.12.21.
//

#include "DEM_Test.hpp"

constexpr int dimension = 2;

typedef DEM_Test<dimension> ParticleMethodType;

int main(int argc, char* argv[]) {

//    openfpm_init(&argc,&argv);

    // Particle data container
    ParticleData<ParticleMethodType> particleData;

    // Transition algorithm
    TransitionCellList<ParticleMethodType> transition(particleData);

    while (!transition.stop(particleData)) {
        transition.run_step(particleData);
    }

//    openfpm_finalize();

    return 0;
}
