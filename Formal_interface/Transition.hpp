//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_TRANSITION_HPP
#define OPENFPM_PDATA_TRANSITION_HPP

#include <Vector/vector_dist.hpp>
#include "ParticleData.hpp"
#include "Particle.hpp"

template <typename ParticleMethodType>
class Transition {

protected:

    typedef typename ParticleMethodType::particleType ParticleType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    ParticleMethodType particleMethod;

//    explicit Transition(ParticleMethod<ParticleType> particleMethod_in) : particleMethod(particleMethod_in), particleData() {}

    int iteration = 0;

    void executeInitialization(ParticleData<ParticleMethodType> &particleData) {
        size_t sz[ParticleMethodType::spaceDimension];
        std::fill(std::begin(sz), std::end(sz), 10);
        auto it2 = particleData.vd.getGridIterator(sz);
        while (it2.isNext())
        {
            particleData.vd.add();
            auto node = it2.get();
            for (int i = 0; i < ParticleMethodType::spaceDimension; i++) {
                particleData.vd.getLastPos()[i] = node.get(i) * it2.getSpacing(i);
                particleData.vd.template getLastProp<0>()[i] = node.get(i);
                particleData.vd.template getLastProp<1>()[i] = 5.0f;
            }
            ++it2;
        }
    }

    void executeEvolution(ParticleData<ParticleMethodType> &particleData) {
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
//            Particle<ParticleType> particle(particleData, p);
            ParticleDataReference<dimension, PositionType, ParticleType> particleDataReference(particleData.vd);
            ParticleRef<dimension, PositionType, ParticleType> particle_ref(particleDataReference, p);

            // call (overriden) evolve method
            particleMethod.evolve(particle_ref);
            ++it2;
        }

//        particleData.vd.template ghost_get<>();

    }

    virtual void executeInteraction(ParticleData<ParticleMethodType> &particleData) {
/*
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particleData, p);

            auto it = particleData.vd.getDomainAndGhostIterator();
            while (it.isNext()) {
                Particle<ParticleType> neighbor(particleData, it.get());
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

    void initialize(ParticleData<ParticleMethodType> &particleData) {
        executeInitialization(particleData);
//        particleData.vd.map();
//        particleData.vd.template ghost_get<0, 1>();
    }

    void run_step(ParticleData<ParticleMethodType> &particleData) {
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

    bool stop(ParticleData<ParticleMethodType> &particleData) {
        return particleMethod.stop();
    }
};

template <typename ParticleMethodType>
class TransitionCellList : public Transition<ParticleMethodType>{
    using typename Transition<ParticleMethodType>::ParticleType;

    CELL_MEMBAL(ParticleMethodType::spaceDimension, float) cellList;

    void executeInteraction(ParticleData<ParticleMethodType> &particleData) override {
/*
        particleData.vd.template updateCellList(cellList);

        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particleData, p);

            auto it = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(particleData.vd.getPos(p)));

//            auto it = this->particleData.vd.getDomainAndGhostIterator();
            while (it.isNext()) {
                Particle<ParticleType> neighbor(particleData, it.get());
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
    explicit TransitionCellList(ParticleData<ParticleMethodType> &particleData) : Transition<ParticleMethodType>(), cellList(particleData.vd.template getCellList<CELL_MEMBAL(ParticleMethodType::spaceDimension, float)>(0.5)) {
        this->initialize(particleData);

    }

};


#endif //OPENFPM_PDATA_TRANSITION_HPP
