//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_DEM_TEST_HPP
#define OPENFPM_PDATA_DEM_TEST_HPP

#include <array>
#include <Vector/vector_dist.hpp>
#include "../Particle.hpp"
#include "../ParticleData.hpp"
#include "../ParticleMethod.hpp"
#include "../Transition.hpp"
#include "../SimulationParameters.hpp"
#include "../InitialCondition.hpp"
#include "../Neighborhood.hpp"


struct DEM_ParticleSignature {
    static constexpr int dimension = 2;
    typedef float position;
    typedef aggregate<float[dimension], float[dimension]> properties;
    typedef FREE_PARTICLES dataStructure;
};

// Property identifier
constexpr int velocity = 0;
constexpr int acceleration = 1;


template <typename ParticleSignature>
class DEM_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    struct GlobalVariable {
        float dt = 0.04;
        float t = 0;
        float t_final = 10.1;
        float r_cut = 0.5;
        float damp = 0.9;
    } globalvar;


    void evolve(Particle<ParticleSignature> particle) override {

        // Apply change of velocity
        particle.template property_test<velocity>() += particle.template property_test<acceleration>();

        // Reset change of velocity
        particle.template property_test<acceleration>() = 0.0f;

        // Euler time-stepping move particles
        Point<dimension, PositionType> delta_velocity = particle.template property<velocity>();
        delta_velocity = delta_velocity * globalvar.dt;
        particle.position_test() += delta_velocity;

/*
        for (int i = 0; i < dimension; i++) {

            // Apply change of velocity
            particle.template property<velocity>()[i] += particle.template property<acceleration>()[i];

            // Reset change of velocity
            particle.template property<acceleration>()[i] = 0.0f;

            // Euler time-stepping move particles
            particle.template position()[i] += particle.template property<velocity>()[i] * globalvar.dt;
        }
*/


    }

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {

        // Declare particle property variables
        Point<dimension, PositionType> p_pos = particle.position();
        Point<dimension, PositionType> n_pos = neighbor.position();
        Point<dimension, PositionType> p_vel = particle.template property<velocity>();
        Point<dimension, PositionType> n_vel = neighbor.template property<velocity>();


        // Check cutoff radius
        if (p_pos.distance(n_pos) > globalvar.r_cut)
            return;

        // Check if particles are moving towards each other
        Point<dimension, PositionType> p_move = p_pos + (p_vel * globalvar.dt);
        Point<dimension, PositionType> n_move = n_pos + (n_vel * globalvar.dt);
        float dist = p_pos.distance2(n_pos);
        float dist_move = (p_move).distance2(n_move);
        if (dist < dist_move)
            return;

        // Compute collision vector
        Point<dimension, PositionType> diff = p_pos - n_pos;
        Point<dimension, PositionType> diff_scaled = diff / p_pos.distance2(n_pos);
        Point<dimension, PositionType> p_collision = (diff * p_vel) * diff_scaled;
        Point<dimension, PositionType> n_collision = (diff * n_vel) * diff_scaled;
        Point<dimension, PositionType> diff_collision = n_collision - p_collision;

        diff_collision = diff_collision * globalvar.damp;

        // Apply collision to particle acceleration
        particle.template property_test<acceleration>() += diff_collision;
        /*
        for (int i = 0; i < dimension; i++) {
            particle.template property<acceleration>()[i] += diff_collision[i];
        }
*/
    }


    bool stop() override {
        std::cout << globalvar.t << std::endl;

        // Check simulation time
        if (globalvar.t > globalvar.t_final)
            return true;

        globalvar.t += globalvar.dt;

        return false;
    }
};

template <typename ParticleSignatureType>
class DEM_SimulationParams : public SimulationParameters<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    // Domain
    Point<dimension, PositionType> domainMin;
    Point<dimension, PositionType> domainMax;

    // Boundary conditions
    periodicity<dimension> boundaryConditions;

//    size_t meshSize[dimension] = {10, 10, 10};

    DEM_SimulationParams() : domainMin(0.0f), domainMax(20.0f) {
        std::fill(std::begin(boundaryConditions.bc), std::end(boundaryConditions.bc), PERIODIC);
    }

/*
    // Mesh initial condition
    typedef InitialConditionMesh initialCondition;
    constexpr static size_t meshSize[dimension] = {18, 18};
*/

    // Random initial condition
    typedef InitialConditionRandom initialCondition;
    int numberParticles = 40;

    // Neighborhood method
    typedef CellListNeighborhood neighborhoodDetermination;
//    float cellWidth = particleMethod.globalvar.r_cut;
    float cellWidth = 0.5;

    void initialization(Particle<ParticleSignatureType> particle) override {

        // Randomize velocity (normal distribution)
        for (int i = 0; i < dimension; i++) {
            particle.template property<velocity>()[i] = this->normalDistribution(0, 3);
        }
    }

};


#endif //OPENFPM_PDATA_DEM_TEST_HPP
