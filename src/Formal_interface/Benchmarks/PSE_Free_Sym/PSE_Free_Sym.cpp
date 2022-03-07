//
// Created by landfried on 02.03.22.
//

#include "PSE_Free_Sym.hpp"

#include <chrono>
#include <fstream>


typedef Benchmark_ParticleMethod1<PSE_ParticleSignature> ParticleMethodType1;
typedef Benchmark_SimulationParams<PSE_ParticleSignature> SimulationParametersType;



int main(int argc, char* argv[]) {


    if (argc >= 2) {
        globalvar.domainSize = std::stof(argv[1]);
        globalvar.meshSize = std::stoi(argv[2]);
        globalvar.recalculate();
    }

    openfpm_init(&argc,&argv);

    auto & vcl = create_vcluster();
    size_t numCPU = vcl.getProcessingUnits();

    std::ofstream output_file;
    std::string filename = "Benchmark_PSE_Free_Sym_" + std::to_string(numCPU) + "core.csv";
    output_file.open(filename);

    output_file << numCPU << " cores,\n";
    output_file << "Time in ms,\n";

    for (int cycle = 0; cycle < 15; ++cycle) {


        // Particle container
        ParticleData<ParticleMethodType1, SimulationParametersType> particleData;

        // State transition
        Transition<ParticleMethodType1, SimulationParametersType> transition(particleData);

        globalvar.t = 0;

        auto t_start = std::chrono::high_resolution_clock::now();

        // Main loop
        while (!transition.stop(particleData)) {
            // Execute simulation step
            transition.run_step(particleData);
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = t_end - t_start;
        std::cout << std::endl << "Duration " << duration.count() << std::endl;

        if (vcl.getProcessUnitID() == 0)
            output_file << (int)duration.count() << ",\n";
    }

    openfpm_finalize();


    output_file.close();

    return 0;
}
