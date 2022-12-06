//
// Created by landfried on 02.03.22.
//
// source /home/peter/openfpm_vars
// make run 

#ifndef OPENFPM_PDATA_SPH_INSTANCE_HPP
#define OPENFPM_PDATA_SPH_INSTANCE_HPP


#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"
#include "Formal_interface/Instance.hpp"
#include "SPH_Algorithm.hpp"

    double GlobalVariable::t = 0;
    double GlobalVariable::dt=0.00016;
    double GlobalVariable::endT= 20.0;

    double GlobalVariable::particleSpacing=1.0/16.0;
    double GlobalVariable::particleSpacingWater=particleSpacing;
    double GlobalVariable::mass=pow(particleSpacing,3)*1000;
    Point<SPH_ParticleSignature::dimension, SPH_ParticleSignature::position>  GlobalVariable::gravity{0.0,0.0,-9.81};
    double GlobalVariable::c0=45.0;
    double GlobalVariable::density0=1000;
    int GlobalVariable::gamma=7;
    double GlobalVariable::nu=1.0/10000.0;
    double GlobalVariable::h=1.3*GlobalVariable::particleSpacing;//characteristic length

    double GlobalVariable::phase=0;
    int GlobalVariable::support=2;
    double GlobalVariable::rc=GlobalVariable::support*GlobalVariable::h;//cutof radius

    double GlobalVariable::domain_min[3] = {-.9, -.9, -.9};
    double GlobalVariable::domain_max[3] = {3.3, 1.3, 1.3};

    // calculate number of particles in each dimension
    size_t GlobalVariable::sz[3] ={
                    uint((GlobalVariable::domain_max[0] - GlobalVariable::domain_min[0]) / GlobalVariable::particleSpacing),
                    uint((GlobalVariable::domain_max[1] - GlobalVariable::domain_min[1]) / GlobalVariable::particleSpacing),
                    uint((GlobalVariable::domain_max[2] - GlobalVariable::domain_min[2]) / GlobalVariable::particleSpacing)};


class DamBreak_SimulationParams : public SimulationParameters<SPH_ParticleSignature> {

public:

    DamBreak_SimulationParams() {

        // Domain
        this->domainMin[0] = g.domain_min[0];
        this->domainMin[1] = g.domain_min[1];
        this->domainMin[2] = g.domain_min[2];
        this->domainMax[0] = g.domain_max[0];
        this->domainMax[1] = g.domain_max[1];
        this->domainMax[2] = g.domain_max[2];

        this->setBoundaryConditions(PERIODIC);
        this->setCutoffRadius(g.rc);
        this->setCellWidth(g.rc);
    }


    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;

    static const int interactionType = INTERACTION_SYMMETRIC;

    // Output
    bool writeOutput = true;
    int writeIteration = 100;

    // DLB
    bool dynamicLoadBalancing = true;

};


class DamBreak_Instance : Instance<SPH_ParticleMethod, DamBreak_SimulationParams> {

    static constexpr int dimension = SPH_ParticleSignature::dimension;
    using PositionType = typename SPH_ParticleSignature::position;


public:

    DamBreak_Instance(ParticleData<SPH_ParticleMethod, DamBreak_SimulationParams> &particleData_in) :
            Instance<SPH_ParticleMethod, DamBreak_SimulationParams>(particleData_in){}

    virtual void shapePlacement() {

        // fluid particles
        PointType waterblockMin{0.15,0.1,0.1};
        PointType waterblockMax{0.55,0.9,0.9};

        auto iterator_fluid = boxIterator(waterblockMin, waterblockMax, g.sz);

        while (iterator_fluid.isNext()) {

            this->addParticle();

            this->position()[0] = iterator_fluid.get().get(0);
            this->position()[1] = iterator_fluid.get().get(1);
            this->position()[2] = iterator_fluid.get().get(2);

            this->property<boundary>() = false;
            this->property<velocity>()[0] = 0.0;
            this->property<velocity>()[1] = 0.0;
            this->property<velocity>()[2] = 0.0;
            this->property<density>() = 1000.0;
            this->property<deltaVelocity>()[0] = 0.0;
            this->property<deltaVelocity>()[1] = 0.0;
            this->property<deltaVelocity>()[2] = 0.0;
            this->property<deltaDensity>() = 0.0;

            ++iterator_fluid;

        }


        // obstacle column

        PointType columnMin{2.0, 3.0/8.0 ,0.0};
        PointType columnMax{2.25,5.0/8.0,1.0};

        auto iterator_column = skinIterator(columnMin, columnMax, g.particleSpacing, g.sz);

        while (iterator_column.isNext()) {

            this->addParticle();

            this->position()[0] = iterator_column.get().get(0);
            this->position()[1] = iterator_column.get().get(1);
            this->position()[2] = iterator_column.get().get(2);

            this->property<boundary>() = true;
            this->property<velocity>()[0] = 0.0;
            this->property<velocity>()[1] = 0.0;
            this->property<velocity>()[2] = 0.0;
            this->property<density>() = 1000.0;
            this->property<deltaVelocity>()[0] = 0.0;
            this->property<deltaVelocity>()[1] = 0.0;
            this->property<deltaVelocity>()[2] = 0.0;
            this->property<deltaDensity>() = 0.0;

            ++iterator_column;

        }


        // pool walls

        PointType poolMin{0.0,0.0,0.0};
        PointType poolMax{3.0,1.0,1.0};

        auto iterator_pool = skinOpenIterator(poolMin, poolMax, g.particleSpacing, g.sz);

        while (iterator_pool.isNext()) {

            this->addParticle();

            this->position()[0] = iterator_pool.get().get(0);
            this->position()[1] = iterator_pool.get().get(1);
            this->position()[2] = iterator_pool.get().get(2);

            this->property<boundary>() = true;
            this->property<velocity>()[0] = 0.0;
            this->property<velocity>()[1] = 0.0;
            this->property<velocity>()[2] = 0.0;
            this->property<density>() = 1000.0;
            this->property<deltaVelocity>()[0] = 0.0;
            this->property<deltaVelocity>()[1] = 0.0;
            this->property<deltaVelocity>()[2] = 0.0;
            this->property<deltaDensity>() = 0.0;

            ++iterator_pool;

        }


    }

    virtual void freePlacement() {}

    void initialization(Particle<SPH_ParticleSignature> particle) override {}


    };



#endif //OPENFPM_PDATA_SPH_INSTANCE_HPP
