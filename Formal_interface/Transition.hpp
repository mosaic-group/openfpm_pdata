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
#include <boost/hana.hpp>

template <typename ParticleMethodType, typename SimulationParametersType>
class Transition {

protected:

//    typedef typename ParticleMethodType::propertyType PropertyType;
//    typedef typename ParticleMethodType::positionType PositionType;
//    static constexpr int dimension = ParticleMethodType::spaceDimension;

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;
    using ParticleDataStructure = typename ParticleSignatureType::dataStructure;


    ParticleMethodType particleMethod;
    SimulationParametersType simulationParameters;

    InitialCondition_Impl<typename SimulationParametersType::initialCondition, ParticleMethodType, SimulationParametersType> initialConditionImplementation;
    Interaction_Impl<typename SimulationParametersType::neighborhoodDetermination, ParticleMethodType, SimulationParametersType> interactionImplementation;

//    explicit Transition(ParticleMethod<PropertyType> particleMethod_in) : particleMethod(particleMethod_in), particleData() {}

    int iteration = 0;



    void executeEvolution(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
//        std::cout << "evolve " << std::endl;
        auto it2 = particleData.getOpenFPMContainer().getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);
//            Particle_VectorDist<dimension, PositionType, PropertyType> particle(particleData.getOpenFPMContainer(), p);

            // call (overriden) evolve method
            particleMethod.evolve(particle);
            ++it2;
        }

    }

    void executeInitialization(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        auto it2 = particleData.getParticleIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleSignatureType> particle(particleData.dataContainer, p);
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
        particleData.getOpenFPMContainer().template map/*<KillParticleWithWarning>*/();

        executeInitialization(particleData);
    }




    void run_step(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
/*
        auto & vcl = create_vcluster();
        if (vcl.getProcessUnitID() == 0) {
            std::cout << "Iteration " << iteration << std::endl;
        }*/

        particleData.getOpenFPMContainer().map();

//        particleData.getOpenFPMContainer().template ghost_get<0, 1>();



        // synchronize ghost for all properties
        particleData.ghost_get_all();



        // call interact method
        interactionImplementation.executeInteraction(particleData);

        // call evolve method
        executeEvolution(particleData);

        // call evolve method (global variable)
        particleMethod.evolveGlobalVariable();


        // write particle data to file
        particleData.getDataContainer().deleteGhost();
        particleData.getDataContainer().write_frame("particles",iteration);

        iteration++;
    }

    bool stop(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        return particleMethod.stop();
    }
};


#endif //OPENFPM_PDATA_TRANSITION_HPP
