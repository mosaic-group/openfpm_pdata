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
float GlobalVariable::t_final = 3;
float GlobalVariable::r_cut = 0.5;
float GlobalVariable::domainSize = 20.0;
bool GlobalVariable::addparticles = true;

class PEC_SimulationParams : public SimulationParameters<DEM_ParticleSignature> {

public:

    PEC_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
    }

    // Random initial condition
    typedef INITIALCONDITION_RANDOM initialCondition;
    int numberParticles = 50;

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;
    float cellWidth = globalvar.r_cut;

    bool writeOutput = true;
    int writeIteration = 30;

};

class PEC_Instance : Instance<DEM_ParticleMethod, PEC_SimulationParams> {

public:

    PEC_Instance(ParticleData<DEM_ParticleMethod, PEC_SimulationParams> &particleData_in) :
    Instance<DEM_ParticleMethod, PEC_SimulationParams>(particleData_in){}

    void initialization(Particle<DEM_ParticleSignature> particle) override {
        // Randomize velocity (normal distribution)

        PARTICLE_WRITE(velocity) = this->normalDistributionPoint(0.0f, 1.0f);

    }

    void freePlacement() {}

    void shapePlacement() {}

};

#endif //OPENFPM_PDATA_EC_INSTANCE_HPP
