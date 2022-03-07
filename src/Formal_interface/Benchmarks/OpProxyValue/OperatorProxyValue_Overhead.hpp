//
// Created by landfried on 08.02.22.
//

#ifndef OPENFPM_PDATA_OPERATORPROXYVALUE_OVERHEAD_HPP
#define OPENFPM_PDATA_OPERATORPROXYVALUE_OVERHEAD_HPP



#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Neighborhood.hpp"


struct Benchmark_ParticleSignature {
    static constexpr int dimension = 2;
    typedef float position;
    typedef aggregate<float[dimension], float[dimension]> properties;
    typedef FREE_PARTICLES dataStructure;
};

// Property identifier
constexpr int prop1 = 0;
constexpr int prop2 = 1;

struct GlobalVariable {
    float domainSize = 20.0;
    int numberCalculations = 15000;
} globalvar;


template <typename ParticleSignature>
class Benchmark_ParticleMethod1 : public ParticleMethod<ParticleSignature> {
    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    void evolve(Particle<ParticleSignature> particle) override {

        particle.template property_test<prop1>() = 0;

        for (int i = 0; i < globalvar.numberCalculations; i++) {
            particle.template property_test<prop1>() += particle.template property_test<prop2>() * i;

        }

    }
};

template <typename ParticleSignature>
class Benchmark_ParticleMethod2 : public ParticleMethod<ParticleSignature> {
    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    void evolve(Particle<ParticleSignature> particle) override {

        for (int d = 0; d < dimension; ++d) {
            particle.template property<prop1>()[d] = 0;
        }

        for (int i = 0; i < globalvar.numberCalculations; i++) {
            for (int d = 0; d < dimension; ++d) {
                particle.template property<prop1>()[d] += particle.template property<prop2>()[d] * i;
            }
        }


    }
};

template <typename ParticleSignatureType>
class Benchmark_SimulationParams : public SimulationParameters<ParticleSignatureType> {
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    Benchmark_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
        this->setNumberParticles(1000);
    }

    // Random initial condition
    typedef INITIALCONDITION_RANDOM initialCondition;

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;

    void initialization(Particle<ParticleSignatureType> particle) override {

        for (int i = 0; i < dimension; i++) {
            particle.template property<prop1>()[i] = 0;
            particle.template property<prop2>()[i] = 1;
        }
    }

};


#endif //OPENFPM_PDATA_OPERATORPROXYVALUE_OVERHEAD_HPP
