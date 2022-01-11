//
// Created by landfried on 21.12.21.
//

#ifndef OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP
#define OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"

template <typename ParticleMethodType>
class SimulationParametersCompiletime {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

public:

//    virtual void initialization(ParticleRef<ParticleMethodType::dimension, typename ParticleMethodType::PositionType, typename ParticleMethodType::ParticleType> particle) {}

    constexpr static PositionType domainMin[dimension] = {0.0, 0.0};
    constexpr static PositionType domainMax[dimension] = {1.0, 1.0};

    // Boundary conditions
    constexpr static size_t boundaryCondition = PERIODIC;

    // Initial condition
    typedef InitialConditionRandom initialCondition;
    constexpr static size_t meshSize[dimension] = {2, 2};
    constexpr static int numberParticles = 1;

};

template <typename ParticleMethodType>
class SimulationParameters {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

public:

//    virtual void initialization(ParticleRef<ParticleMethodType::dimension, typename ParticleMethodType::PositionType, typename ParticleMethodType::ParticleType> particle) {}

    PositionType domainMin[dimension] = {0.0, 0.0};
    PositionType domainMax[dimension] = {1.0, 1.0};

    // Boundary conditions
    size_t boundaryCondition = PERIODIC;

    // Initial condition
    typedef InitialConditionRandom initialCondition;
    size_t meshSize[dimension] = {2, 2};
    int numberParticles = 1;

};

#endif //OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP
