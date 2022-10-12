//
// Created by landfried on 17.12.21.
//

#include "FI_parallel_version.hpp"
#include <chrono>
#include <fstream>

typedef PSE_ParticleMethod<PSE_ParticleSignature> ParticleMethodType;
typedef PSE_SimulationParams<PSE_ParticleSignature> SimulationParametersType;

int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);


    auto & vcl = create_vcluster();
    size_t numCPU = vcl.getProcessingUnits();

    std::ofstream output_file;
    std::string filename = "Benchmark_FI_parallel_" + std::to_string(numCPU) + "core.csv";
    output_file.open(filename);

    output_file << numCPU << " cores,\n";
    output_file << "Time in ms,\n";

    for (int cycle = 0; cycle < 15; ++cycle) {

        // Particle container
        ParticleData<ParticleMethodType, SimulationParametersType> particleData;

        // State transition
        Transition<ParticleMethodType, SimulationParametersType> transition(particleData);

        globalvar.t = 0;

        auto t_start2 = std::chrono::high_resolution_clock::now();

        // Main loop
        while (!transition.stop(particleData)) {

            // Execute simulation step
            transition.run_step(particleData);
        }

        auto t_end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration2 = t_end2 - t_start2;
        std::cout << std::endl << "Duration " << duration2.count() << std::endl;

        if (vcl.getProcessUnitID() == 0)
            output_file << (int)duration2.count() << ",\n";


    }

    output_file.close();


    openfpm_finalize();

    return 0;
}