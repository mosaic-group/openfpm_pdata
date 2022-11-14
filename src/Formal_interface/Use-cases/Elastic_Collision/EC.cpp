//
// Created by landfried on 27.09.22.
//

#include "EC_Algorithm.hpp"
#include "EC_Instance.hpp"

typedef EC_ParticleMethod<EC_ParticleSignature> ParticleMethodType;
typedef EC_SimulationParams<EC_ParticleSignature> SimulationParametersType;
typedef EC_Instance<ParticleMethodType, SimulationParametersType> InstanceType;

int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    // Particle container
    ParticleData<ParticleMethodType, SimulationParametersType> particleData;

    // State transition
    Transition<ParticleMethodType, SimulationParametersType, InstanceType> transition(particleData);

    // Main loop
    while (!transition.stop(particleData)) {

        // Execute simulation step
        transition.run_step(particleData);
    }

    openfpm_finalize();

    return 0;
}
