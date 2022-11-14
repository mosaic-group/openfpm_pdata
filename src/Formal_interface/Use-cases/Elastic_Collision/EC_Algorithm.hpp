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
#include "Formal_interface/Util.hpp"


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

        // Compute collision
        auto diff = neighbor.position() - particle.position();
        PARTICLE(acceleration) += diff / abs2(diff) * scalarProduct(diff, NEIGHBOR(velocity) - PARTICLE(velocity));
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
