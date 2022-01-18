//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_TRANSITION_HPP
#define OPENFPM_PDATA_TRANSITION_HPP

#include <Vector/vector_dist.hpp>
#include "ParticleData.hpp"
#include "Particle.hpp"
#include "InitialCondition.hpp"
#include "Neighborhood.hpp"
#include <random>

template <typename ParticleMethodType, typename SimulationParametersType>
class Transition {

protected:

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    InitialCondition_Impl<typename SimulationParametersType::initialCondition, ParticleMethodType, SimulationParametersType> initialConditionImplementation;
    Interaction_Impl<typename SimulationParametersType::neighborhoodDetermination, ParticleMethodType, SimulationParametersType> interactionImplementation;

//    explicit Transition(ParticleMethod<PropertyType> particleMethod_in) : particleMethod(particleMethod_in), particleData() {}

    int iteration = 0;



    void executeEvolution(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
//        std::cout << "evolve " << std::endl;
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<dimension, PositionType, PropertyType> particle(particleData.vd, p);

            // call (overriden) evolve method
            particleMethod.evolve(particle);
            ++it2;
        }

    }

    void executeInitialization(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<dimension, PositionType, PropertyType> particle(particleData.vd, p);
            simulationParameters.initialization(particle);
            ++it2;
        }

    }

public:

    Transition(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) : interactionImplementation(particleData) {
        initializeParticles(particleData);
    }

    void initializeParticles(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        initialConditionImplementation.initialization(particleData);
        particleData.vd.template map<KillParticleWithWarning>();

        executeInitialization(particleData);
    }

    void run_step(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
/*
        auto & vcl = create_vcluster();
        if (vcl.getProcessUnitID() == 0) {
            std::cout << "Iteration " << iteration << std::endl;
        }*/

        particleData.vd.map();

        particleData.vd.template ghost_get<0, 1>();

//        executeInteraction(particleData);
        interactionImplementation.executeInteraction(particleData);
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


#endif //OPENFPM_PDATA_TRANSITION_HPP
