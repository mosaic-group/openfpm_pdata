// Created by Johannes Pahkle
//

#include "SPH_Instance.hpp"
#include "SPH_Algorithm.hpp"
#include <chrono>

int main(int argc, char* argv[]) {

    openfpm_init(&argc, &argv);

    // Particle container
    ParticleData<SPH_ParticleMethod, DamBreak_SimulationParams> particleData;

    // State transition
    Transition<SPH_ParticleMethod, DamBreak_SimulationParams, DamBreak_Instance> transition(particleData);

    // Main loop
    while (!transition.stop(particleData)) {

        // Execute simulation step
        transition.run_step(particleData);
    }

    openfpm_finalize();

    return 0;
}
