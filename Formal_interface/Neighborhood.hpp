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

// Primary template
template <typename NeighborhoodDetermination, typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl {
public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}
};



template <typename ParticleMethodType, typename SimulationParametersType>
class Interaction_Impl<AllParticlesNeighborhood, ParticleMethodType, SimulationParametersType> {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    ParticleMethodType particleMethod;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {}


    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        // iterate through all particles
        auto iteratorAll = particleData.vd.getDomainIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<dimension, PositionType, PropertyType> particle(particleData.vd, p);

            // iterate through all particles as neighbors
            auto iteratorNeighbors = particleData.vd.getDomainIterator();
            while (iteratorNeighbors.isNext()) {
                auto n = iteratorNeighbors.get();
                Particle<dimension, PositionType, PropertyType> neighbor(particleData.vd, n);

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

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    CELL_MEMBAL(dimension, float) cellList;

public:

    explicit Interaction_Impl(ParticleData<ParticleMethodType, SimulationParametersType> &particleData): cellList(particleData.vd.template getCellList<CELL_MEMBAL(ParticleMethodType::spaceDimension, float)>(simulationParameters.cellWidth)) {}

    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        particleData.vd.template updateCellList(cellList);

        // iterate through all particles
        auto iteratorAll = particleData.vd.getDomainIterator();
        while (iteratorAll.isNext())
        {
            auto p = iteratorAll.get();
            Particle<dimension, PositionType, PropertyType> particle(particleData.vd, p);

            // iterate through all particles in neighbor cells
            auto iteratorNeighbors = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(particleData.vd.getPos(p)));

            while (iteratorNeighbors.isNext()) {
                auto n = iteratorNeighbors.get();
                Particle<dimension, PositionType, PropertyType> neighbor(particleData.vd, n);

                if (particle != neighbor) {
                    particleMethod.interact(particle, neighbor);
                }
                ++iteratorNeighbors;
            }
            ++iteratorAll;
        }
    }

};


#endif //OPENFPM_PDATA_INTERACTION_HPP
