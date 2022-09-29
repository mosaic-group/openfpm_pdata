//
// Created by landfried on 27.09.22.
//

#ifndef OPENFPM_PDATA_EC_INSTANCE_HPP
#define OPENFPM_PDATA_EC_INSTANCE_HPP

#include <array>
#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"
#include "Formal_interface/Use-cases/Elastic_Collision/EC_Algorithm.hpp"


float GlobalVariable::dt = 0.001;
float GlobalVariable::t = 0;
float GlobalVariable::t_final = 50;
float GlobalVariable::r_cut = 0.5;
float GlobalVariable::domainSize = 20.0;

template <typename ParticleSignatureType>
class EC_SimulationParams : public SimulationParameters<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    // Domain
    Point<dimension, PositionType> domainMin;
    Point<dimension, PositionType> domainMax;

    EC_SimulationParams() : domainMin(0.0f), domainMax(globalvar.domainSize) {
        this->setBoundaryConditions(PERIODIC);
    }

    // Random initial condition
    typedef INITIALCONDITION_RANDOM initialCondition;
    int numberParticles = 50;

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;
    float cellWidth = globalvar.r_cut;


    void initialization(Particle<ParticleSignatureType> particle) override {

        // Randomize velocity (normal distribution)
        for (int i = 0; i < dimension; i++) {
            particle.template property<velocity>()[i] = this->normalDistribution(0, 1);
        }
    }

    bool writeOutput = true;
    int writeIteration = 25;

};


#endif //OPENFPM_PDATA_EC_INSTANCE_HPP
