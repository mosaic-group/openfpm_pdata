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
double GlobalVariable::t_final = 100;

double GlobalVariable::domainSize = 40.0;
int GlobalVariable::meshSize = 128;
double GlobalVariable::meshSpacing = GlobalVariable::domainSize / GlobalVariable::meshSize;
double GlobalVariable::epsilon = GlobalVariable::meshSpacing;
double GlobalVariable::r_cut = 3 * GlobalVariable::epsilon;
double GlobalVariable::D = 0.01;
//double GlobalVariable::kernelFactor = GlobalVariable::D * 15.0 * pow(GlobalVariable::meshSpacing, 3) / (pow(GlobalVariable::epsilon, 5)  * pow(M_PI, 2));


class Diffusion_SimulationParams : public SimulationParameters<PSE_ParticleSignature> {


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

    int writeIteration = 100;

};


class Diffusion_Instance : Instance<PSE_ParticleMethod, Diffusion_SimulationParams> {

    static constexpr int dimension = PSE_ParticleSignature::dimension;

public:

    Diffusion_Instance(ParticleData<PSE_ParticleMethod, Diffusion_SimulationParams> &particleData_in) :
            Instance<PSE_ParticleMethod, Diffusion_SimulationParams>(particleData_in){}

    void initialization(Particle<PSE_ParticleSignature> particle) override {

        PointType centerPoint = globalvar.domainSize / 4;
        if (particle.position() == centerPoint)
            PARTICLE(concentration) = 1 / pow(globalvar.meshSpacing, 3);

    }

    void freePlacement() {}

    void shapePlacement() {}


};


#endif //OPENFPM_PDATA_DIF_INSTANCE_HPP
