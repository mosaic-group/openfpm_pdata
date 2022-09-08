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
#include "Formal_interface/Neighborhood.hpp"
#include "Formal_interface/Instance.hpp"
#include "SPH_Algorithm.hpp"



template <typename ParticleSignatureType>
class SPH_SimulationParams : public SimulationParameters<ParticleSignatureType> {
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    SPH_SimulationParams() {

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

    // DLB
    bool dynamicLoadBalancing = true;



};


class Instance1 : Instance<SPH_ParticleMethod<SPH_ParticleSignature>, SPH_SimulationParams<SPH_ParticleSignature>> {

    static constexpr int dimension = SPH_ParticleSignature::dimension;
    using PositionType = typename SPH_ParticleSignature::position;


public:

    Instance1(ParticleData<SPH_ParticleMethod<SPH_ParticleSignature>, SPH_SimulationParams<SPH_ParticleSignature>> &particleData_in) :
            Instance<SPH_ParticleMethod<SPH_ParticleSignature>, SPH_SimulationParams<SPH_ParticleSignature>>(particleData_in){}

    virtual void shapePlacement() {

        // fluid particles
        Point<dimension, PositionType> waterblockMin{0.15,0.1,0.1};
        Point<dimension, PositionType> waterblockMax{0.55,0.9,0.9};

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

        Point<dimension, PositionType> columnMin{2.0, 3.0/8.0 ,0.0};
        Point<dimension, PositionType> columnMax{2.25,5.0/8.0,1.0};

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

        Point<dimension, PositionType> poolMin{0.0,0.0,0.0};
        Point<dimension, PositionType> poolMax{3.0,1.0,1.0};

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

    virtual void freePlacement() {


        // double particleSpacing=1.0/64.0;
        // double particleSpacingWater=particleSpacing;

        //     GlobalVariable g;
        //     g.endT=20;
        //     g.t=0.0;
        //     g.dt=0.00005;

        //     g.mass=power(particleSpacing,3)*1000;
        //     g.h=1.3*particleSpacing;

        //     g.epsilon=0.01;
        //     g.gamma=7;
        //     g.gravity=Point<dimension, PositionType>(0.0,0.0,-9.81);
        //     g.density0=1000;
        //     g.c0=45.0;
        //     g.nu=1.0/10000.0;
        //     g.phase=0;

        //     int support=2;
        //     g.rc=support*g.h;


        //these should be a multible of the particleSpacing
        Point<dimension, PositionType> poolMin{0.0,0.0,0.0};
        Point<dimension, PositionType> poolMax{3.0,1.0,1.0};
        Point<dimension, PositionType> columnMin{2.0, 3.0/8.0 ,0.0};
        Point<dimension, PositionType> columnMax{2.25,5.0/8.0,1.0};
        Point<dimension, PositionType> waterblockMin{0.0,0.0,0.0};
        Point<dimension, PositionType> waterblockMax{0.5,1.0,1.0};
        // waterblockMin+=Point<dimension, PositionType>(g.particleSpacing);
        // waterblockMax-=Point<dimension, PositionType>(g.particleSpacing);

        // g.domain_min=poolMin-Point<dimension, PositionType>(g.particleSpacing*g.support);
        // g.domain_max=Point<dimension, PositionType>{poolMax[0]+g.particleSpacing*g.support,
        //                                 poolMax[1]+g.particleSpacing*g.support,
        //                                 poolMax[2]
        // };



/*

        // fluid
        this->addParticle();
        this->position()[0] = 0.0;
        this->position()[1] = 0.0;
        this->position()[2] = 1.0;
        this->property<boundary>() = false;
        this->property<velocity>()[0] = 0.0;
        this->property<velocity>()[1] = 0.0;
        this->property<velocity>()[2] = 0.0;
        this->property<density>() = 1000.0;
        this->property<deltaVelocity>()[0] = 0.0;
        this->property<deltaVelocity>()[1] = 0.0;
        this->property<deltaVelocity>()[2] = 0.0;
        this->property<deltaDensity>() = 0.0;

        this->addParticle();
        this->position()[0] = 0.0;
        this->position()[1] = 0.0;
        this->position()[2] = 2.0;
        this->property<boundary>() = false;
        this->property<velocity>()[0] = 0.0;
        this->property<velocity>()[1] = 0.0;
        this->property<velocity>()[2] = 0.0;
        this->property<density>() = 1000.0;
        this->property<deltaVelocity>()[0] = 0.0;
        this->property<deltaVelocity>()[1] = 0.0;
        this->property<deltaVelocity>()[2] = 0.0;
        this->property<deltaDensity>() = 0.0;

        // wall
        this->addParticle();
        this->position()[0] = 0.0;
        this->position()[1] = 0.0;
        this->position()[2] = 0.0;
        this->property<boundary>() = true;
        this->property<velocity>()[0] = 0.0;
        this->property<velocity>()[1] = 0.0;
        this->property<velocity>()[2] = 0.0;
        this->property<density>() = 1000.0;
        this->property<deltaVelocity>()[0] = 0.0;
        this->property<deltaVelocity>()[1] = 0.0;
        this->property<deltaVelocity>()[2] = 0.0;
        this->property<deltaDensity>() = 0.0;
*/




        Point<dimension, PositionType> wallMin;
        Point<dimension, PositionType> wallMax;

        double d= g.particleSpacingWater/2.0;

/*

// fluid block
        wallMin[0]=waterblockMin[0]+g.particleSpacing;
        wallMin[1]=waterblockMin[1]+g.particleSpacing;
        wallMin[2]=waterblockMin[2]+g.particleSpacing;
        wallMax[0]=waterblockMax[0]-g.particleSpacing;
        wallMax[1]=waterblockMax[1]-g.particleSpacing;
        wallMax[2]=waterblockMax[2]-g.particleSpacing;
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = false;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

*/





/*

        //Pool Boundary
        d=g.particleSpacing/2.0;

        //Botom Pool Boundary 
        wallMin[0]=poolMin[0]-g.particleSpacing*g.support;
        wallMin[1]=poolMin[1]-g.particleSpacing*g.support;
        wallMin[2]=poolMin[2]-g.particleSpacing*g.support;
        wallMax[0]=poolMax[0]+g.particleSpacing*g.support;
        wallMax[1]=poolMax[1]+g.particleSpacing*g.support;
        wallMax[2]=poolMin[2]-g.particleSpacing;
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Left Pool Boundary 
        wallMin[0]=poolMin[0]-g.particleSpacing*g.support;
        wallMin[1]=poolMin[1]-g.particleSpacing*g.support;
        wallMin[2]=poolMin[2];
        wallMax[0]=poolMin[0]-g.particleSpacing;
        wallMax[1]=poolMax[1]+g.particleSpacing*g.support;
        wallMax[2]=poolMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Right Pool Boundary 
        wallMin[0]=poolMax[0]+g.particleSpacing;
        wallMin[1]=poolMin[1]-g.particleSpacing*g.support;
        wallMin[2]=poolMin[2];
        wallMax[0]=poolMax[0]+g.particleSpacing*g.support;
        wallMax[1]=poolMax[1]+g.particleSpacing*g.support;
        wallMax[2]=poolMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Front Pool Boundary 
        wallMin[0]=poolMin[0];
        wallMin[1]=poolMin[1]-g.particleSpacing*g.support;
        wallMin[2]=poolMin[2];
        wallMax[0]=poolMax[0];
        wallMax[1]=poolMin[1]-g.particleSpacing;
        wallMax[2]=poolMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Back Pool Boundary 
        wallMin[0]=poolMin[0];
        wallMin[1]=poolMax[1]+g.particleSpacing;
        wallMin[2]=poolMin[2];
        wallMax[0]=poolMax[0];
        wallMax[1]=poolMax[1]+g.particleSpacing*g.support;
        wallMax[2]=poolMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }


*/








/*        //Top Column Boundary
        wallMin[0]=columnMin[0]-g.particleSpacing*g.support;
        wallMin[1]=columnMin[1]-g.particleSpacing*g.support;
        wallMin[2]=columnMax[2]+g.particleSpacing;
        wallMax[0]=columnMax[0]+g.particleSpacing*g.support;
        wallMax[1]=columnMax[1]+g.particleSpacing*g.support;
        wallMax[2]=columnMax[2]+g.particleSpacing*g.support;
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Left Column Boundary 
        wallMin[0]=columnMin[0]-g.particleSpacing*g.support;
        wallMin[1]=columnMin[1]-g.particleSpacing*g.support;
        wallMin[2]=columnMin[2];
        wallMax[0]=columnMin[0]-g.particleSpacing;
        wallMax[1]=columnMax[1]+g.particleSpacing*g.support;
        wallMax[2]=columnMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Right Column Boundary 
        wallMin[0]=columnMax[0]+g.particleSpacing;
        wallMin[1]=columnMin[1]-g.particleSpacing*g.support;
        wallMin[2]=columnMin[2];
        wallMax[0]=columnMax[0]+g.particleSpacing*g.support;
        wallMax[1]=columnMax[1]+g.particleSpacing*g.support;
        wallMax[2]=columnMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Front Column Boundary 
        wallMin[0]=columnMin[0];
        wallMin[1]=columnMin[1]-g.particleSpacing*g.support;
        wallMin[2]=columnMin[2];
        wallMax[0]=columnMax[0];
        wallMax[1]=columnMin[1]-g.particleSpacing;
        wallMax[2]=columnMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }

        //Back Column Boundary 
        wallMin[0]=columnMin[0];
        wallMin[1]=columnMax[1]+g.particleSpacing;
        wallMin[2]=columnMin[2];
        wallMax[0]=columnMax[0];
        wallMax[1]=columnMax[1]+g.particleSpacing*g.support;
        wallMax[2]=columnMax[2];
        for (double x = wallMin[0]; x <= wallMax[0]+d; x+=g.particleSpacing) {
            for (double y = wallMin[1]; y <= wallMax[1]+d; y+=g.particleSpacing) {
                for (double z = wallMin[2]; z <= wallMax[2]+d; z+=g.particleSpacing) {
                    this->addParticle();
                    this->position()[0] = x;
                    this->position()[1] = y;
                    this->position()[2] = z;
                    this->property<boundary>() = true;
                    this->property<velocity>()[0] = 0.0;
                    this->property<velocity>()[1] = 0.0;
                    this->property<velocity>()[2] = 0.0;
                    this->property<density>() = 1000.0;
                    this->property<deltaVelocity>()[0] = 0.0;
                    this->property<deltaVelocity>()[1] = 0.0;
                    this->property<deltaVelocity>()[2] = 0.0;
                    this->property<deltaDensity>() = 0.0;
                }
            }
        }*/
    }

/*
    virtual void setGlobalVariable() {

        g.t=0;
        g.dt=0.00004;
        g.endT= 2.0;

        g.particleSpacing=1.0/64.0;
        g.particleSpacingWater = g.particleSpacing;
        g.mass=pow(g.particleSpacing,3)*1000;
        g.gravity[0] = 0.0;
        g.gravity[1] = 0.0;
        g.gravity[2] = -9.81;
        g.c0=45.0;
        g.density0=1000;
        g.gamma=7;
        g.nu=1.0/10000.0;
        g.h=1.3*g.particleSpacing;//characteristic length

        g.phase=0;
        g.support=2;
        g.rc=g.support*g.h;//cutof radius
        g.epsilon=0.01;

        g.domain_min[0] = -.9;
        g.domain_min[1] = -.9;
        g.domain_min[2] = -.9;

        g.domain_max[0] = 3.3;
        g.domain_max[1] = 1.3;
        g.domain_max[2] = 1.3;


        // calculate number of particles in each dimension
        g.sz[0] = uint((g.domain_max[0] - g.domain_min[0]) / g.particleSpacing);
        g.sz[1] = uint((g.domain_max[1] - g.domain_min[1]) / g.particleSpacing);
        g.sz[2] = uint((g.domain_max[2] - g.domain_min[2]) / g.particleSpacing);


    }
*/

};



#endif //OPENFPM_PDATA_SPH_INSTANCE_HPP
