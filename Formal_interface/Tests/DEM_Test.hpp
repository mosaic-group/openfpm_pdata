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


struct GlobalVariable {
    float dt = 0.04;
    float t = 0;
    float t_final = 30.1;
    float r_cut = 0.5;
    float damp = 0.9;
    float domainSize = 20.0;
} globalvar;


template <typename ParticleSignature>
class DEM_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:


    void evolve(Particle<ParticleSignature> particle) override {

        // Apply change of velocity
        particle.template property_test<velocity>() += particle.template property_test<acceleration>();

        // Reset change of velocity
        particle.template property_test<acceleration>() = 0.0f;

        // Euler time-stepping move particles
        particle.position_test() += particle.template property_test<velocity>() * globalvar.dt;

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

    }

    void evolveGlobalVariable() override {

        // advance time
        globalvar.t += globalvar.dt;
        std::cout << "\r" << int(globalvar.t / globalvar.t_final * 100) << "%" << std::flush;

    }

    bool stop() override {

        // Check simulation time
        if (globalvar.t > globalvar.t_final)
            return true;

        return false;
    }
};

template <typename ParticleSignatureType>
class DEM_SimulationParams : public SimulationParameters<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

/*
    size_t domainMin[dimension] = {-1, 0};
    size_t domainMax[dimension] = {3, 5};
*/

    DEM_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
        this->setCutoffRadius(globalvar.r_cut);
        this->setNumberParticles(40);
        this->setCellWidth(globalvar.r_cut);
    }



/*
    // Mesh initial condition
    typedef INITIALCONDITION_MESH initialCondition;
    constexpr static size_t meshSize[dimension] = {18, 18};
*/

    // Random initial condition
    typedef INITIALCONDITION_RANDOM initialCondition;

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;
//    typedef NEIGHBORHOOD_ALLPARTICLES neighborhoodDetermination;

    void initialization(Particle<ParticleSignatureType> particle) override {

        // Randomize velocity (normal distribution)
        for (int i = 0; i < dimension; i++) {
            particle.template property<velocity>()[i] = this->normalDistribution(0, 3);
        }
    }

};


#endif //OPENFPM_PDATA_DEM_TEST_HPP
