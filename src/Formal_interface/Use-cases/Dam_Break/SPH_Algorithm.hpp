//
// Created by landfried on 02.03.22.
//
// source /home/peter/openfpm_vars
// make run 

#ifndef OPENFPM_PDATA_SPH_ALGORITHM_HPP
#define OPENFPM_PDATA_SPH_ALGORITHM_HPP


#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"
#include "Formal_interface/Instance.hpp"
#include "Formal_interface/Util.hpp"
#include "Formal_interface/Alias.hpp"

// Property identifier
constexpr int positionOld = 0;
constexpr int velocity = 1;
constexpr int velocityOld = 2;
constexpr int density = 3;
constexpr int densityOld = 4;
constexpr int deltaVelocity = 5;
constexpr int deltaDensity = 6;
constexpr int boundary = 7;


struct SPH_ParticleSignature {
    static constexpr int dimension = 3;
    typedef double position;
    typedef aggregate<double[dimension], double[dimension], double[dimension], double, double, double[dimension],double, bool> properties;
    typedef FREE_PARTICLES dataStructure;
};


struct GlobalVariable {

    static double t;
    static double dt;
    static double endT;

    static double particleSpacing;
    static double particleSpacingWater;
    static double mass;
    static Point<SPH_ParticleSignature::dimension, SPH_ParticleSignature::position> gravity;
    static double c0;
    static double density0;
    static int gamma;
    static double nu;
    static double h;//characteristic length

    static double phase;
    static int support;
    static double rc;//cutof radius

    static double domain_min[3];
    static double domain_max[3];

    // calculate number of particles in each dimension
    static size_t sz[3];

} g;


class SPH_ParticleMethod : public ParticleMethod<SPH_ParticleSignature> {
    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;


    double pressure_density2(double density){
        return g.c0*g.c0*g.density0/g.gamma * (pow(density/g.density0, g.gamma)-1)/density/density;
    }

public:


    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {

        PointType r_pq = neighbor.position() - particle.position();
        double f_pq=  pow(1.0 - distance(neighbor, particle) / 2.0 / g.h, 3);
        double vr = scalarProduct(r_pq, NEIGHBOR(velocity) - PARTICLE(velocity));

        // Compute change of velocity
        double interim01 = pressure_density2(NEIGHBOR(density)) + pressure_density2(PARTICLE(density));
        double interim02 = 10.0 * g.nu / distance2(neighbor, particle) * vr;
        PARTICLE(deltaVelocity) += (interim01 - interim02 / PARTICLE(density)) * r_pq * f_pq;
        NEIGHBOR(deltaVelocity) -= (interim01 - interim02 / NEIGHBOR(density)) * r_pq * f_pq;

        // Compute change of density
        double densityChange= vr * f_pq;
        PARTICLE(deltaDensity) += densityChange;
        NEIGHBOR(deltaDensity) += densityChange;
    }

    void evolve(Particle<ParticleSignature> particle) override {
        double prefact = g.mass* -5.0*21.0/16.0/M_PI/pow(g.h,5);

        PointType acceleration = g.gravity + PARTICLE(deltaVelocity) * prefact;

        double densityAcceleration = PARTICLE(deltaDensity) * prefact;

        if (g.phase==0){
            if (PARTICLE(boundary) == false){
                // fluid

                // move particles half step
                PARTICLE(positionOld) = particle.position();
                particle.position() += g.dt/2.0f * PARTICLE(velocity);

                // change velocity half step
                PARTICLE(velocityOld) = PARTICLE(velocity);
                PARTICLE(velocity) += g.dt/2.0 * acceleration;
            }

            // fluid + boundary

            // change density
            PARTICLE(densityOld) = PARTICLE(density);
            PARTICLE(density) += g.dt/2.0 * densityAcceleration;
        }
        else {
            if (PARTICLE(boundary) == false) {
                // fluid

                // move particle from original position
                PointType acceleration_halfstep = g.dt  / 2.0 * acceleration;
                PointType velocity_corrected = g.dt * (PARTICLE(velocityOld) + acceleration_halfstep);
                particle.position() = PARTICLE(positionOld) + velocity_corrected;

                // change velocity
                PointType acceleration_corrected = g.dt * acceleration;
                PARTICLE(velocity) = PARTICLE(velocityOld) + acceleration_corrected;
            }

            // fluid + boundary

            // change density
            PARTICLE(density) = PARTICLE(densityOld) + g.dt * densityAcceleration;
        }

        //set to 0 to have a fresh accumulators for the next time step
        PARTICLE(deltaVelocity) = PointType (0.0);
        PARTICLE(deltaDensity) = 0.0;
    }

    void evolveGlobalVariable() override {

        if (g.phase==0){
            g.phase=1;
        }
        else{
            g.t+=g.dt;//advance the current time by the time step
            g.phase=0;
        }

    }

    bool stop() override {

        // Check simulation time
        if (g.t > g.endT)
            return true;

        return false;
    }
};




#endif //OPENFPM_PDATA_SPH_ALGORITHM_HPP
