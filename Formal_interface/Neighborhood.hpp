//
// Created by landfried on 13.01.22.
//

#ifndef OPENFPM_PDATA_INTERACTION_HPP
#define OPENFPM_PDATA_INTERACTION_HPP

#include <iostream>
#include "ParticleData.hpp"
#include "Particle.hpp"

struct AllParticlesNeighborhood {};
struct CellListNeighborhood {};
struct MeshNeighborhood {};

// Primary template
template <typename NeighborhoodDetermination, typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl {
public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}
};



template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<AllParticlesNeighborhood, ParticleMethodType, SimulationParametersType> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}


    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        // iterate through all particles
        auto iteratorAll = particleData.getContainer().getDomainIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles as neighbors
            auto iteratorNeighbors = particleData.getContainer().getDomainIterator();
            while (iteratorNeighbors.isNext()) {
                auto n = iteratorNeighbors.get();
                Particle<ParticleSignatureType> neighbor(particleData.dataContainer, n);

                if (particle != neighbor) {
                    particleMethod.interact(particle, neighbor);
                }
                ++iteratorNeighbors;
            }
            ++iteratorAll;
        }
    }
};


template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<CellListNeighborhood, ParticleMethodType, SimulationParametersType> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    CELL_MEMBAL(dimension, float) cellList;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellList(particleData.getContainer().template getCellList<CELL_MEMBAL(dimension, float)>(simulationParameters.cellWidth)) {}

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        particleData.getContainer().template updateCellList(cellList);

        // iterate through all particles
        auto iteratorAll = particleData.getContainer().getDomainIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles in neighbor cells
            auto iteratorNeighbors = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(particleData.getContainer().getPos(p)));

            while (iteratorNeighbors.isNext()) {
                auto n = iteratorNeighbors.get();
                Particle<ParticleSignatureType> neighbor(particleData.dataContainer, n);

                if (particle != neighbor) {
                    particleMethod.interact(particle, neighbor);
                }
                ++iteratorNeighbors;
            }
            ++iteratorAll;
        }
    }

};

template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<MeshNeighborhood, ParticleMethodType, SimulationParametersType> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}


    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        // iterate through all particles
        auto iteratorAll = particleData.getContainer().getDomainIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all neighbor mesh nodes
            for (int d = 0; d < dimension; ++d) {
                auto nMinus = p.move(d, -1);
                Particle<ParticleSignatureType> neighborMinus(particleData.dataContainer, nMinus);
                particleMethod.interact(particle, neighborMinus);

                auto nPlus = p.move(d, 1);
                Particle<ParticleSignatureType> neighborPlus(particleData.dataContainer, nPlus);
                particleMethod.interact(particle, neighborPlus);
            }

            ++iteratorAll;
        }
    }
};


#endif //OPENFPM_PDATA_INTERACTION_HPP
