//
// Created by landfried on 27.09.22.
//

#ifndef OPENFPM_PDATA_DIF_ALGORITHM_HPP
#define OPENFPM_PDATA_DIF_ALGORITHM_HPP


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

struct PSE_ParticleSignature {
    static constexpr int dimension = 3;
    typedef double position;
    typedef aggregate<size_t, double, double> properties;
    typedef MESH_PARTICLES dataStructure;
};

// Property identifier
constexpr int concentration = 1;
constexpr int accumulator = 2;


struct GlobalVariable {
    static double dt;
    static double t;
    static double t_final;

    static double domainSize;
    static int meshSize;
    static double meshSpacing;
    static double epsilon;
    static double r_cut;
    static double D;
    static double kernel;
} globalvar;

template <typename ParticleSignature>
class PSE_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    void evolve(Particle<ParticleSignature> particle) override {
        PARTICLE(concentration) += PARTICLE(accumulator) * globalvar.kernel;
        PARTICLE(accumulator) = 0;
    }

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {
        PositionType r_pq2 = distance2(particle, neighbor);
        double exchange = (NEIGHBOR(concentration) - PARTICLE(concentration))
                          / (1 + pow(r_pq2 / globalvar.epsilon / globalvar.epsilon, 5)) ;
        PARTICLE(accumulator) += exchange;
        NEIGHBOR(accumulator) -= exchange;
    }

    void evolveGlobalVariable() override {

        // advance time
        globalvar.t += globalvar.dt;
    }

    bool stop() override {

        // Check simulation time
        if (globalvar.t > globalvar.t_final)
            return true;

        return false;
    }
};


#endif //OPENFPM_PDATA_DIF_ALGORITHM_HPP
