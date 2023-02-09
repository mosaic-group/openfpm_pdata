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
#include <random>

float GlobalVariable::dt = 0.001;
float GlobalVariable::t = 0;
float GlobalVariable::t_final = 30.002;
float GlobalVariable::r_cut = 0.70;
float GlobalVariable::domainSize = 20.0;
//bool GlobalVariable::addparticles = true;

class PEC_SimulationParams : public SimulationParameters<DEM_ParticleSignature> {

public:

    PEC_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
        //setMeshSize(10);
    }

    // Random initial condition
    // typedef INITIALCONDITION_RANDOM initialCondition;
    typedef INITIALCONDITION_RANDOM initialCondition;

    float cutoff_radius = globalvar.r_cut;

    int numberParticles = 1000;

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;
    float cellWidth = globalvar.r_cut;

    bool writeOutput = true;
    int writeIteration = 10;

};

class PEC_Instance : Instance<DEM_ParticleMethod, PEC_SimulationParams> {

public:

    PEC_Instance(ParticleData<DEM_ParticleMethod, PEC_SimulationParams> &particleData_in) :
    Instance<DEM_ParticleMethod, PEC_SimulationParams>(particleData_in){}

    void initialization(Particle<DEM_ParticleSignature> particle) override {
        // Randomize velocity (normal distribution)

        PARTICLE(velocity) = 
        // {(float)rand() / (float)RAND_MAX - 0.5,
        //                      (float)rand() / (float)RAND_MAX - 0.5,
        //                      (float)rand() / (float)RAND_MAX - 0.5};
        this->normalDistributionPoint(0.0f, 1.0f);
        
        PARTICLE(deltaVelocity) = {0.0,0.0};
    }

    void freePlacement() {}

    void shapePlacement() {}

};

#endif //OPENFPM_PDATA_EC_INSTANCE_HPP
