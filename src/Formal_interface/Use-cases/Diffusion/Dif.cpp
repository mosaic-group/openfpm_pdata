//
// Created by landfried on 27.09.22.
//

#include "Dif_Algorithm.hpp"
#include "Dif_Instance.hpp"

typedef PSE_ParticleMethod<PSE_ParticleSignature> ParticleMethodType;
typedef Diffusion_SimulationParams<PSE_ParticleSignature> SimulationParametersType;
typedef Diffusion_Instance<ParticleMethodType, SimulationParametersType> InstanceType;

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
