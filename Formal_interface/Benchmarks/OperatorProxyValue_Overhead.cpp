//
// Created by landfried on 08.02.22.
//

#include "OperatorProxyValue_Overhead.hpp"
#include <chrono>
#include <fstream>


typedef Benchmark_ParticleMethod1<Benchmark_ParticleSignature> ParticleMethodType1;
typedef Benchmark_ParticleMethod2<Benchmark_ParticleSignature> ParticleMethodType2;
typedef Benchmark_SimulationParams<Benchmark_ParticleSignature> SimulationParametersType;


constexpr int TIMESTEPS = 5000;

int main(int argc, char* argv[]) {


    std::ofstream output_file;
    output_file.open("Benchmark_OpProxyValue.csv");

    output_file << "Time in ms\n";
    output_file << "Forloop,OperationProxyValue,\n";

    openfpm_init(&argc,&argv);

    auto & vcl = create_vcluster();

    for (int cycle = 0; cycle < 15; ++cycle) {

        // ---- using for loop ----

        // Particle container
        ParticleData<ParticleMethodType2, SimulationParametersType> particleData2;

        // State transition
        Transition<ParticleMethodType2, SimulationParametersType> transition2(particleData2);

        auto t_start2 = std::chrono::high_resolution_clock::now();

        // Main loop
        for (int i = 0; i < TIMESTEPS; ++i) {
            // Execute simulation step
            transition2.run_step(particleData2);
        }

        auto t_end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration2 = t_end2 - t_start2;
        std::cout << std::endl << "Duration (forloop) " << duration2.count() << std::endl;

        if (vcl.getProcessUnitID() == 0)
            output_file << duration2.count() << ",";



        // ---- using OperationProxyValue ----

        // Particle container
        ParticleData<ParticleMethodType1, SimulationParametersType> particleData;

        // State transition
        Transition<ParticleMethodType1, SimulationParametersType> transition(particleData);

        auto t_start = std::chrono::high_resolution_clock::now();

        // Main loop
        for (int i = 0; i < TIMESTEPS; ++i) {
            // Execute simulation step
            transition.run_step(particleData);
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = t_end - t_start;
        std::cout << std::endl << "Duration (opproxy) " << duration.count() << std::endl;

        if (vcl.getProcessUnitID() == 0)
            output_file << duration.count() << ",\n";



    }

    openfpm_finalize();


    output_file.close();

    return 0;
}
