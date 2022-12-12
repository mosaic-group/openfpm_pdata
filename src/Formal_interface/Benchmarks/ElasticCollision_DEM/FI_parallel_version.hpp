//
// Created by landfried on 17.12.21.
//

#ifndef OPENFPM_PDATA_PSE_TEST_HPP
#define OPENFPM_PDATA_PSE_TEST_HPP


#include <array>
#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"


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
    float dt = 0.001;
    float t = 0;
    float t_final = .5;
    float r_cut = 0.08;
    float damp = 0.9;
    float domainSize = 9.0;
    int number_particles = 200000;
} globalvar;


template <typename ParticleSignature>
class DEM_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:


    void evolve(Particle<ParticleSignature> particle) override {

        // Apply change of velocity
        particle.template property<velocity>() += particle.template property<acceleration>();

        // Reset change of velocity
        particle.template property<acceleration>() = 0.0f;

        // Euler time-stepping move particles
        particle.position() += particle.template property<velocity>() * globalvar.dt;

    }

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {

        // Declare particle property variables
        Point<dimension, PositionType> p_pos = particle.position_raw();
        Point<dimension, PositionType> n_pos = neighbor.position_raw();
        Point<dimension, PositionType> p_vel = particle.template property_raw<velocity>();
        Point<dimension, PositionType> n_vel = neighbor.template property_raw<velocity>();


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
        Point<dimension, PositionType> diff = n_pos - p_pos;
        Point<dimension, PositionType> diff_scaled = diff / n_pos.distance2(p_pos);
        Point<dimension, PositionType> diff_vel = n_vel - p_vel;

//        Point<dimension, PositionType> diff_collision = diff_scaled * scalarProduct(diff_vel, diff);
        Point<dimension, PositionType> diff_collision = diff * scalarProduct(diff_scaled, diff_vel);

        diff_collision = diff_collision; //* globalvar.damp;


        // Apply collision to particle acceleration
        particle.template property<acceleration>() += diff_collision;


    }


    bool stop() override {
//        std::cout << "\r" << int(globalvar.t / globalvar.t_final * 100) << "%" << std::flush;

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


//    size_t meshSize[dimension] = {10, 10, 10};

    DEM_SimulationParams() : domainMin(0.0f), domainMax(globalvar.domainSize) {
        this->setBoundaryConditions(PERIODIC);
    }

/*
    // Mesh initial condition
    typedef INITIALCONDITION_MESH initialCondition;
    constexpr static size_t meshSize[dimension] = {18, 18};
*/

    // Random initial condition
    typedef INITIALCONDITION_RANDOM initialCondition;
    int numberParticles = globalvar.number_particles;

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;
//    float cellWidth = particleMethod.globalvar.r_cut;
    float cellWidth = globalvar.r_cut;


    void initialization(Particle<ParticleSignatureType> particle) override {

        // Randomize velocity (normal distribution)
        for (int i = 0; i < dimension; i++) {
            particle.template property<velocity>()[i] = this->normalDistribution(0, 1);
        }
    }

    bool writeOutput = false;
//    int writeIteration = 25;

};



#endif //OPENFPM_PDATA_PSE_TEST_HPP
