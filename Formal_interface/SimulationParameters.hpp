//
// Created by landfried on 21.12.21.
//

#ifndef OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP
#define OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"


template <typename ParticleSignatureType>
class SimulationParameters {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;
    using ParticleDataStructure = typename ParticleSignatureType::dataStructure;

public:

    SimulationParameters() : gen(rd()) {}


//    virtual void initialization(ParticleRef<ParticleMethodType::dimension, typename ParticleMethodType::PositionType, typename ParticleMethodType::PropertyType> particle) {}

    PositionType domainMin[dimension] = {0.0, 0.0};
    PositionType domainMax[dimension] = {1.0, 1.0};

    // Boundary conditions
    size_t boundaryConditions[dimension] = {PERIODIC, PERIODIC};

    // Initial condition
    typedef InitialConditionRandom initialCondition;
    size_t meshSize[dimension] = {5, 5};
    int numberParticles = 1;

    virtual void initialization(Particle<ParticleSignatureType> particle) {}


    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()

    float normalDistribution(float mean, float stdev) {
        std::normal_distribution<> dis(mean, stdev);
        return dis(gen);
    }
};

#endif //OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP
