//
// Created by landfried on 21.12.21.
//

#ifndef OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP
#define OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"
#include "Constants.hpp"


template <typename ParticleSignatureType>
class SimulationParameters {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;
    using ParticleDataStructure = typename ParticleSignatureType::dataStructure;

public:

    SimulationParameters() : gen(rd()) {}


    // Domain
    PositionType domainMin[dimension] = {0.0};
    PositionType domainMax[dimension] = {1.0};

    // Boundary conditions
    size_t boundaryConditions[dimension] = {PERIODIC};

    // Initial condition
    typedef INITIALCONDITION_RANDOM initialCondition;
    int numberParticles = 1;

    // Mesh
    size_t meshSize[dimension] = {5};
    PositionType meshSpacing = 0.1;

    // Cell list
    PositionType cellWidth = 0.3;

    // Interaction
    typedef NEIGHBORHOOD_ALLPARTICLES neighborhoodDetermination;
    static const int interactionType = INTERACTION_PULL;
    PositionType cutoff_radius = 0.3;



protected:

    // Initialization methods

    void setDomain(PositionType min, PositionType max) {
        std::fill(std::begin(domainMin), std::end(domainMin), min);
        std::fill(std::begin(domainMax), std::end(domainMax), max);
    }

    void setDomain(PositionType max) {
        setDomain(0, max);
    }

    void setBoundaryConditions(size_t value) {
        std::fill(std::begin(boundaryConditions), std::end(boundaryConditions), value);
    }

    void setMeshSize(int value) {
        std::fill(std::begin(meshSize), std::end(meshSize), value);
    }

    void setCutoffRadius(PositionType value) {
        cutoff_radius = value;
    }

    void setMeshSpacing(PositionType value) {
        meshSpacing = value;
    }

    void setNumberParticles(int value) {
        numberParticles = value;
    }

    void setCellWidth(PositionType value) {
        cellWidth = value;
    }

    // Utility methods

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()

    float normalDistribution(float mean, float stdev) {
        std::normal_distribution<> dis(mean, stdev);
//        std::uniform_real_distribution<> dis(0, 10);
        return dis(gen);
    }


public:

    // Particle methods

    virtual void initialization(Particle<ParticleSignatureType> particle) {}

};

#endif //OPENFPM_PDATA_SIMULATIONPARAMETERS_HPP
