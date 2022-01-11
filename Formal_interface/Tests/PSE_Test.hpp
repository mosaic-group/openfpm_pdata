//
// Created by landfried on 17.12.21.
//

#ifndef OPENFPM_PDATA_PSE_TEST_HPP
#define OPENFPM_PDATA_PSE_TEST_HPP

#include <Vector/vector_dist.hpp>
#include "../Particle.hpp"
#include "../ParticleData.hpp"
#include "../ParticleMethod.hpp"
#include "../Transition.hpp"

// Particle properties
static constexpr int concentration = 0;
static constexpr int accumulator = 1;

// Global variable properties
static constexpr int dt = 0;
static constexpr int t = 1;
static constexpr int t_final = 2;
static constexpr int r_cut = 3;
static constexpr int epsilon = 4;
static constexpr int D = 5;


// Position type
using position_type = float;

// Particle type
template <int dimension>
using property_type_n = aggregate<float[dimension], float[dimension]>;

// Global variable type
using globalvar_type = aggregate<float, float, float, float, float, float>;


template <int dimension>
class PSE_Test : ParticleMethod<dimension, position_type, property_type_n<dimension>, globalvar_type> {
public:

    bool meshParticles;
    constexpr static position_type domainMin = 0.0;
    constexpr static position_type domainMax = 1.0;
    constexpr static position_type meshSpacing = 0.005;
//    constexpr static position_type meshDim = 1.0;


    bool stop(GlobalVar<globalvar_type> globalVar) override {
        return globalVar.property<t>() >= globalVar.property<t_final>();
    }

    void interact(GlobalVar<globalvar_type> globalVar, Particle<property_type_n<dimension>> particle, Particle<property_type_n<dimension>> neighbor) override {
        position_type d2 = getDistance2(particle, neighbor);
        particle.property<accumulator>() += (particle.property<concentration>() - neighbor.property<concentration>())
                / (power(d2 / globalVar.property<epsilon>() / globalVar.property<epsilon>(), 5) + 1.0);
    }

    void evolve(GlobalVar<globalvar_type> globalVar, Particle<property_type_n<dimension>> particle, Particle<property_type_n<dimension>> neighbor) override {
        particle.property<concentration>() += globalVar.property<dt>() * globalVar.property<D>() * 15.0 * power(meshSpacing / globalVar.property<epsilon>(), 3)
                / power(globalVar.property<epsilon>() * M_PI, 2) * particle.property<accumulator>();
        particle.property<accumulator>() = 0;
    }

    void evolveGlobalVar(GlobalVar<globalvar_type> globalVar) override {
        globalVar.property<t>() += globalVar.property<dt>();
    }


};


#endif //OPENFPM_PDATA_PSE_TEST_HPP
