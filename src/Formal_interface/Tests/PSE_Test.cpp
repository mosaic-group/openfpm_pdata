//
// Created by landfried on 17.12.21.
//

#include "PSE_Test.hpp"
#include <chrono>
#include <fstream>

typedef PSE_ParticleMethod<PSE_ParticleSignature> ParticleMethodType;
typedef PSE_SimulationParams<PSE_ParticleSignature> SimulationParametersType;

int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    // Particle container
    ParticleData<ParticleMethodType, SimulationParametersType> particleData;

    // State transition
    Transition<ParticleMethodType, SimulationParametersType> transition(particleData);

    auto t_start2 = std::chrono::high_resolution_clock::now();

    // Main loop
    while (!transition.stop(particleData)) {

        // Execute simulation step
        transition.run_step(particleData);
    }

    auto t_end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration2 = t_end2 - t_start2;
    std::cout << std::endl << "Duration " << duration2.count() << std::endl;

    std::ofstream output_file;
    output_file.open("DEM_Sym.csv");
    output_file << duration2.count() << ",\n";
    output_file.close();


    openfpm_finalize();

    return 0;
}