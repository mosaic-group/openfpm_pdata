//
// Created by landfried on 23.02.22.
//

#ifndef OPENFPM_PDATA_CELLLIST_HPP
#define OPENFPM_PDATA_CELLLIST_HPP

#include "ParticleData.hpp"
#include "Particle.hpp"
#include "Constants.hpp"



template <typename ParticleMethodType, typename SimulationParametersType, int InteractionType = SimulationParametersType::interactionType>
class CellList_Impl {};

template <typename ParticleMethodType, typename SimulationParametersType>
class CellList_Impl<ParticleMethodType, SimulationParametersType, INTERACTION_SYMMETRIC> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    using CellListType = typename vector_dist<dimension, PositionType, PropertyType>::CellList_type;

    CellListType cellList;

public:

    explicit CellList_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellList(
            createCellList(particleData)) {

    }

    CellListType createCellList(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        return particleData.getOpenFPMContainer().getCellListSym(simulationParameters.cellWidth);
    }

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        // update symmetric cell list
        particleData.getOpenFPMContainer().updateCellListSym(cellList);

        // iterate through all particles
        auto iteratorAll = particleData.getParticleIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles in neighbor cells
            auto iteratorNeighbors = cellList.template getNNIteratorSym<NO_CHECK>(cellList.getCell(
                    particleData.getOpenFPMContainer().getPos(p)), p.getKey(), particleData.getOpenFPMContainer().getPosVector());

            while (iteratorNeighbors.isNext()) {
                auto n = iteratorNeighbors.get();
                Particle<ParticleSignatureType> neighbor(particleData.dataContainer, n);

                // calculate distance
                Point<dimension, PositionType> n_pos = neighbor.position();
                Point<dimension, PositionType> p_pos = particle.position();
                PositionType distance = p_pos.distance(n_pos);

                // check if distance <= cutoff radius
                if (distance <= simulationParameters.cutoff_radius) {
                    if (particle != neighbor) {
                        particleMethod.interact(particle, neighbor);
                    }
                }

                ++iteratorNeighbors;
            }
            ++iteratorAll;
        }
    }

};

template <typename ParticleMethodType, typename SimulationParametersType>
class CellList_Impl<ParticleMethodType, SimulationParametersType, INTERACTION_PULL> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    using CellListType = typename vector_dist<dimension, PositionType, PropertyType>::CellList_type;

    CellListType cellList;

public:

    explicit CellList_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellList(
            createCellList(particleData)) {

    }

    CellListType createCellList(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        return particleData.getOpenFPMContainer().getCellList(simulationParameters.cellWidth);
    }

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        // update symmetric cell list
        particleData.getOpenFPMContainer().updateCellList(cellList);

        // iterate through all particles
        auto iteratorAll = particleData.getParticleIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles in neighbor cells
            auto iteratorNeighbors = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(
                    particleData.getOpenFPMContainer().getPos(p)));

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
class CellList_Impl<ParticleMethodType, SimulationParametersType, INTERACTION_PUSH> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    using CellListType = typename vector_dist<dimension, PositionType, PropertyType>::CellList_type;

    CellListType cellList;

public:

    explicit CellList_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellList(
            createCellList(particleData)) {

    }

    CellListType createCellList(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        return particleData.getOpenFPMContainer().getCellList(simulationParameters.cellWidth);
    }

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        // update symmetric cell list
        particleData.getOpenFPMContainer().updateCellList(cellList);

        // iterate through all particles
        auto iteratorAll = particleData.getParticleIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles in neighbor cells
            auto iteratorNeighbors = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(
                    particleData.getOpenFPMContainer().getPos(p)));

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
class CellList_Impl<ParticleMethodType, SimulationParametersType, INTERACTION_PULL_PUSH> {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    using CellListType = typename vector_dist<dimension, PositionType, PropertyType>::CellList_type;

    CellListType cellList;

public:

    explicit CellList_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellList(
            createCellList(particleData)) {

    }

    CellListType createCellList(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        return particleData.getOpenFPMContainer().getCellList(simulationParameters.cellWidth);
    }

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        // update symmetric cell list
        particleData.getOpenFPMContainer().updateCellList(cellList);

        // iterate through all particles
        auto iteratorAll = particleData.getParticleIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);

            // iterate through all particles in neighbor cells
            auto iteratorNeighbors = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(
                    particleData.getOpenFPMContainer().getPos(p)));

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


#endif //OPENFPM_PDATA_CELLLIST_HPP
