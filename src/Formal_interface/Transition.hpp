//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_TRANSITION_HPP
#define OPENFPM_PDATA_TRANSITION_HPP

#include "Vector/vector_dist.hpp"
#include "ParticleData.hpp"
#include "Particle.hpp"
#include "InitialCondition.hpp"
#include "Interaction_Impl.hpp"
#include "Instance.hpp"
#include <random>
#include <boost/hana.hpp>

template <typename ParticleMethodType, typename SimulationParametersType, typename InstanceType = Instance<ParticleMethodType, SimulationParametersType>>
class Transition {

protected:

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

    // Dynamic load balancing
    ModelSquare dlb_model;
    int steps_rebalance = simulationParameters.balanceIteration;



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

        InstanceType instance(particleData);


        auto & vcl = create_vcluster();
        // initialize particles on single core
        if (vcl.getProcessUnitID() == 0) {
            instance.freePlacement();
        }


        instance.shapePlacement();

        // place particles
        // random or on a mesh
        initialConditionImplementation.initialization(particleData);

        // particle-wise initialization
        executeInitialization(particleData);


        // dynamic load balancing
        if (simulationParameters.dynamicLoadBalancing) {
            particleData.getOpenFPMContainer().map();
            particleData.getOpenFPMContainer().addComputationCosts(dlb_model);
            particleData.getOpenFPMContainer().getDecomposition().decompose();
        }

        // distribute particles across cores
        particleData.getOpenFPMContainer().map();

        // synchronize ghost for all properties
        particleData.ghost_get_all();

    }


    void run_step(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {

        if (simulationParameters.dynamicLoadBalancing) {
            if (iteration % steps_rebalance == 0) {
                // dynamic load balancing
                particleData.getOpenFPMContainer().map();
                particleData.getOpenFPMContainer().addComputationCosts(dlb_model);
                particleData.getOpenFPMContainer().getDecomposition().redecompose(steps_rebalance);
            }
        }

        // distribute particles across cores
        particleData.getOpenFPMContainer().map();

        // synchronize ghost for all properties
        particleData.ghost_get_all();

        // call interact method
        interactionImplementation.executeInteraction(particleData);

        // synchronize ghost for all properties
//        particleData.ghost_get_all();

        // call evolve method
        executeEvolution(particleData);

        // call evolve method (global variable)
        particleMethod.evolveGlobalVariable();


        // write particle data to file
        if (simulationParameters.writeOutput) {
            if (iteration % simulationParameters.writeIteration == 0)
            {
                // write particles
                particleData.getDataContainer().deleteGhost();
                particleData.getDataContainer().write_frame("particles", iteration);
            }
        }



        iteration++;
    }

    bool stop(ParticleData<ParticleMethodType, SimulationParametersType> &particleData) {
        return particleMethod.stop();
    }
};


#endif //OPENFPM_PDATA_TRANSITION_HPP
