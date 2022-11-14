//
// Created by landfried on 27.09.22.
//

#ifndef OPENFPM_PDATA_EC_ALGORITHM_HPP
#define OPENFPM_PDATA_EC_ALGORITHM_HPP


#include <array>
#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"
#include "Formal_interface/Alias.hpp"


struct EC_ParticleSignature {
    static constexpr int dimension = 2;
    typedef float position;
    typedef aggregate<float[dimension], float[dimension]> properties;
    typedef FREE_PARTICLES dataStructure;
};

// Property identifier
constexpr int velocity = 0;
constexpr int acceleration = 1;


struct GlobalVariable {
    static float dt;
    static float t;
    static float t_final;
    static float r_cut;
    static float domainSize;
} globalvar;

template <typename ParticleSignature>
class EC_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    void evolve(Particle<ParticleSignature> particle) override {

        // Apply change of velocity
        PARTICLE(velocity) += PARTICLE(acceleration);

        // Reset change of velocity
        PARTICLE(acceleration) = 0.0f;

        // Euler time-stepping move particles
        particle.position() += PARTICLE(velocity) * globalvar.dt;

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

        // Compute collision vector
        Point<dimension, PositionType> diff = n_pos - p_pos;
        Point<dimension, PositionType> diff_scaled = diff / n_pos.distance2(p_pos);
        Point<dimension, PositionType> diff_vel = n_vel - p_vel;
        Point<dimension, PositionType> diff_collision = diff_scaled * scalarProduct(diff, diff_vel);

        // Apply collision to particle acceleration
        PARTICLE(acceleration) += diff_collision;
    }


    void evolveGlobalVariable() {
        // advance time
        globalvar.t += globalvar.dt;
    }

    bool stop() override {
        // Check simulation time
        return globalvar.t > globalvar.t_final;
    }
};


#endif //OPENFPM_PDATA_EC_ALGORITHM_HPP
