//
// Created by landfried on 27.09.22.
//

#ifndef OPENFPM_PDATA_DIF_INSTANCE_HPP
#define OPENFPM_PDATA_DIF_INSTANCE_HPP

#include <array>
#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"
#include "Formal_interface/Instance.hpp"
#include "Formal_interface/Use-cases/Diffusion/Dif_Algorithm.hpp"


double GlobalVariable::dt = 0.05;
double GlobalVariable::t = 0;
double GlobalVariable::t_final = 1;

double GlobalVariable::domainSize = 40.0;
int GlobalVariable::meshSize = 64;
double GlobalVariable::meshSpacing = GlobalVariable::domainSize / GlobalVariable::meshSize;
double GlobalVariable::epsilon = GlobalVariable::meshSpacing;
double GlobalVariable::r_cut = 3 * GlobalVariable::epsilon;
double GlobalVariable::D = 0.01;
double GlobalVariable::kernel = GlobalVariable::dt * GlobalVariable::D * 15.0 *
        pow(GlobalVariable::meshSpacing/GlobalVariable::epsilon, 3)  / pow(GlobalVariable::epsilon * M_PI, 2);


template <typename ParticleSignatureType>
class Diffusion_SimulationParams : public SimulationParameters<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    Diffusion_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
        this->setMeshSize(globalvar.meshSize);
        this->setCutoffRadius(globalvar.r_cut);
        this->setMeshSpacing(globalvar.meshSpacing);
    }

    typedef NEIGHBORHOOD_MESH neighborhoodDetermination;
    static const int interactionType = INTERACTION_SYMMETRIC;

    int writeIteration = 10;

};


template <typename ParticleMethodType, typename SimulationParametersType>
class Diffusion_Instance : Instance<ParticleMethodType, SimulationParametersType> {

    static constexpr int dimension = PSE_ParticleSignature::dimension;
    using PositionType = typename PSE_ParticleSignature::position;


public:

    Diffusion_Instance(ParticleData<ParticleMethodType, SimulationParametersType> &particleData_in) :
            Instance<ParticleMethodType, SimulationParametersType>(particleData_in){}

    void initialization(Particle<PSE_ParticleSignature> particle) override {
        bool centerParticle = true;
        for (int i = 0; i < dimension; ++i) {
            if (particle.position()[i] != globalvar.domainSize / 4)
                centerParticle = false;
        }

        if (centerParticle) {
            PARTICLE(concentration) = 1 / pow(globalvar.meshSpacing, 3);
        }
    }

    void freePlacement() {}

    void shapePlacement() {}


};


#endif //OPENFPM_PDATA_DIF_INSTANCE_HPP
