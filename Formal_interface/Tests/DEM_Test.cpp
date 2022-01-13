//
// Created by landfried on 06.12.21.
//

#include "DEM_Test.hpp"

constexpr int dimension = 2;

typedef DEM_ParticleMethod<dimension> ParticleMethodType;
typedef DEM_SimulationParams<ParticleMethodType> SimulationParametersType;

int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    // Particle container
    ParticleData<ParticleMethodType, SimulationParametersType> particleData;

    // State transition
    Transition<ParticleMethodType, SimulationParametersType> transition(particleData);

    // Main loop
    while (!transition.stop(particleData)) {

        // Execute simulation step
        transition.run_step(particleData);
    }

    openfpm_finalize();

    return 0;
}
