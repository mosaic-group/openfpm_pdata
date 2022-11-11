// Created by Johannes Pahkle
//

#include "SPH_Instance.hpp"
#include "SPH_Algorithm.hpp"
#include <chrono>


typedef SPH_ParticleMethod<SPH_ParticleSignature> ParticleMethodType;
typedef DamBreak_SimulationParams<SPH_ParticleSignature> SimulationParametersType;
typedef DamBreak_Instance InstanceType;

int main(int argc, char* argv[]) {

    openfpm_init(&argc, &argv);

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
