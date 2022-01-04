//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_TRANSITION_HPP
#define OPENFPM_PDATA_TRANSITION_HPP

#include <Vector/vector_dist.hpp>
#include "ParticleData.hpp"
#include "Particle.hpp"
#include <random>

template <typename ParticleMethodType, typename InitialConditionType>
class Transition {

protected:

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    ParticleMethodType particleMethod;

//    explicit Transition(ParticleMethod<PropertyType> particleMethod_in) : particleMethod(particleMethod_in), particleData() {}

    int iteration = 0;

    void executeInitialization(ParticleData<ParticleMethodType, InitialConditionType> &particleData) {

        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis_pos(ParticleMethodType::domainMin, ParticleMethodType::domainMax);
        std::uniform_real_distribution<> dis_vel(-0.6, 0.6);

        size_t sz[ParticleMethodType::spaceDimension];
        std::fill(std::begin(sz), std::end(sz), 5);
//        sz = {10, 10, 10};
        auto it2 = particleData.vd.getGridIterator(sz);
        while (it2.isNext())
        {
            particleData.vd.add();
            auto node = it2.get();
            for (int i = 0; i < ParticleMethodType::spaceDimension; i++) {
//                particleData.vd.getLastPos()[i] = node.get(i) * (it2.getSpacing(i)) + dis_pos(gen);
                particleData.vd.getLastPos()[i] = node.get(i) * it2.getSpacing(i);
//                particleData.vd.template getLastProp<0>()[i] = node.get(i);
                particleData.vd.template getLastProp<1>()[i] = dis_vel(gen);
            }
            ++it2;
        }
    }

    void executeEvolution(ParticleData<ParticleMethodType, InitialConditionType> &particleData) {
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

    virtual void executeInteraction(ParticleData<ParticleMethodType, InitialConditionType> &particleData) {
/*
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<PropertyType> particle(particleData, p);

            auto it = particleData.vd.getDomainAndGhostIterator();
            while (it.isNext()) {
                Particle<PropertyType> neighbor(particleData, it.get());
                if (particle != neighbor) {
//                    std::cout << particle.template property<0>() << " neighbor prop 0 " << neighbor.getParticleData().vd.template getProp<0>(neighbor.getID()) << std::endl;
                    particleMethod.interact(particle, neighbor);
                }
                ++it;
            }
            ++it2;
        }
*/
    }


public:

    void initialize(ParticleData<ParticleMethodType, InitialConditionType> &particleData) {
        executeInitialization(particleData);
//        particleData.vd.map();
//        particleData.vd.template ghost_get<0, 1>();
    }

    void run_step(ParticleData<ParticleMethodType, InitialConditionType> &particleData) {
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

    bool stop(ParticleData<ParticleMethodType, InitialConditionType> &particleData) {
        return particleMethod.stop();
    }
};

template <typename ParticleMethodType, typename InitialConditionType>
class TransitionCellList : public Transition<ParticleMethodType, InitialConditionType>{
    using typename Transition<ParticleMethodType, InitialConditionType>::PropertyType;

    CELL_MEMBAL(ParticleMethodType::spaceDimension, float) cellList;

    void executeInteraction(ParticleData<ParticleMethodType, InitialConditionType> &particleData) override {
/*
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
*/
    }

public:
    explicit TransitionCellList(ParticleData<ParticleMethodType, InitialConditionType> &particleData) : Transition<ParticleMethodType, InitialConditionType>(), cellList(particleData.vd.template getCellList<CELL_MEMBAL(ParticleMethodType::spaceDimension, float)>(0.5)) {
        this->initialize(particleData);

    }

};


#endif //OPENFPM_PDATA_TRANSITION_HPP
