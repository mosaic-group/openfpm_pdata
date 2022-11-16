//
// Created by landfried on 27.09.22.
//

#include "Dif_Algorithm.hpp"
#include "Dif_Instance.hpp"

int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    // Particle container
    ParticleData<PSE_ParticleMethod, Diffusion_SimulationParams> particleData;

    // State transition
    Transition<PSE_ParticleMethod, Diffusion_SimulationParams, Diffusion_Instance> transition(particleData);

    // Main loop
    while (!transition.stop(particleData)) {

        // Execute simulation step
        transition.run_step(particleData);
    }

    openfpm_finalize();

    return 0;
}
