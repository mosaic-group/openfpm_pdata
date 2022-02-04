//
// Created by landfried on 17.12.21.
//

#include "PSE_Test.hpp"

typedef PSE_ParticleMethod<PSE_ParticleSignature> ParticleMethodType;
typedef PSE_SimulationParams<PSE_ParticleSignature> SimulationParametersType;

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