//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_TRANSITION_HPP
#define OPENFPM_PDATA_TRANSITION_HPP

#include <Vector/vector_dist.hpp>
#include "ParticleData.hpp"
#include "Particle.hpp"
#include "InitialCondition.hpp"
#include <random>

template <typename ParticleMethodType, typename SimulationParametersType>
class Transition {

protected:

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    ParticleMethodType particleMethod;
    InitialCondition_Impl<typename SimulationParametersType::initialCondition, ParticleMethodType, SimulationParametersType> initialConditionImplementation;

//    explicit Transition(ParticleMethod<PropertyType> particleMethod_in) : particleMethod(particleMethod_in), particleData() {}

    int iteration = 0;

    virtual void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

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
//                    std::cout << particle.template property<0>() << " neighbor prop 0 " << neighbor.getParticleData().vd.template getProp<0>(neighbor.getID()) << std::endl;
//                    std::cout << "CellList" << std::endl;
                    particleMethod.interact(particle, neighbor);
                }
                ++iteratorNeighbors;
            }
            ++iteratorAll;
        }
    }

    void executeEvolution(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<dimension, PositionType, PropertyType> particle(particleData.vd, p);

            // call (overriden) evolve method
            particleMethod.evolve(particle);
            ++it2;
        }

//        particleData.vd.template ghost_get<>();

    }


public:

    void initialize(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        initialConditionImplementation.initialization(particleData);
        particleData.vd.map();

    }

    void run_step(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
/*
        auto & vcl = create_vcluster();
        if (vcl.getProcessUnitID() == 0) {
            std::cout << "Iteration " << iteration << std::endl;
        }*/

        particleData.vd.map();

        particleData.vd.template ghost_get<0, 1>();

        executeInteraction(particleData);
        executeEvolution(particleData);

        particleData.vd.deleteGhost();
        particleData.vd.write_frame("particles",iteration);

//        std::cout << "iteration " << iteration << std::endl;
        iteration++;
    }

    bool stop(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        return particleMethod.stop();
    }
};

template <typename ParticleMethodType, typename SimulationParametersType>
class TransitionCellList : public Transition<ParticleMethodType, SimulationParametersType>{
    using typename Transition<ParticleMethodType, SimulationParametersType>::PropertyType;

    CELL_MEMBAL(ParticleMethodType::spaceDimension, float) cellList;

/*
    void executeInteraction(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) override {
        particleData.vd.template updateCellList(cellList);

        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<PropertyType> particle(particleData, p);

            auto it = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(particleData.vd.getPos(p)));

//            auto it = this->particleData.vd.getDomainAndGhostIterator();
            while (it.isNext()) {
                Particle<PropertyType> neighbor(particleData, it.get());
                if (particle != neighbor) {
//                    std::cout << particle.template property<0>() << " neighbor prop 0 " << neighbor.getParticleData().vd.template getProp<0>(neighbor.getID()) << std::endl;
//                    std::cout << "CellList" << std::endl;
                    this->particleMethod.interact(particle, neighbor);
                }
                ++it;
            }
            ++it2;
        }
    }
*/

public:
    explicit TransitionCellList(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) : Transition<ParticleMethodType, SimulationParametersType>(), cellList(particleData.vd.template getCellList<CELL_MEMBAL(ParticleMethodType::spaceDimension, float)>(0.5)) {
        this->initialize(particleData);

    }

};


#endif //OPENFPM_PDATA_TRANSITION_HPP
