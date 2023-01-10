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
    typedef aggregate<double, double> properties;
    typedef MESH_PARTICLES dataStructure;
};

// Property identifier
constexpr int concentration = 0;
constexpr int accumulator = 1;


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
} globalvar;

class PSE_ParticleMethod : public ParticleMethod<PSE_ParticleSignature> {

public:

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {
        double exchange = (NEIGHBOR(concentration) - PARTICLE(concentration)) /
                (1.0 + pow(distance2(neighbor, particle), 5));
        PARTICLE(accumulator) += exchange;
        NEIGHBOR(accumulator) -= exchange;
    }

    void evolve(Particle<ParticleSignature> particle) override {

        PARTICLE(concentration) += globalvar.dt * globalvar.D * 15.0 *
                pow(globalvar.meshSpacing, 3) / pow(globalvar.epsilon * M_PI, 2) * PARTICLE(accumulator);
        PARTICLE(accumulator) = 0;
    }

    void evolveGlobalVariable() override {

        // advance time
        globalvar.t += globalvar.dt;
    }

    bool stop() override {

        // Check simulation time
        return globalvar.t > globalvar.t_final;
    }
};


#endif //OPENFPM_PDATA_DIF_ALGORITHM_HPP
