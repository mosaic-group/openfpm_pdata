//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_DEM_TEST_HPP
#define OPENFPM_PDATA_DEM_TEST_HPP

#include <Vector/vector_dist.hpp>
#include "../Particle.hpp"
#include "../ParticleData.hpp"
#include "../ParticleMethod.hpp"
#include "../Transition.hpp"
#include "../GlobalVar.hpp"


// Position type
using position_type = float;

// n-dimensional properties
template <int dimension>
using property_type_n = aggregate<float[dimension], float[dimension]>;

typedef aggregate<float, float, float, float> globalvar_type;

template <int dimension>
class DEM_Test : public ParticleMethod<dimension, position_type, property_type_n<dimension>> {
    //, globalvar_type
public:
    static constexpr int position = 0;
    static constexpr int velocity = 1;

    static constexpr int dt = 0;
    static constexpr int time = 1;
    static constexpr int t_final = 2;
    static constexpr int r_cut = 3;


    void evolve(GlobalVar<globalvar_type> globalVar, Particle<dimension, position_type, property_type_n<dimension>> particle) override {

        // Euler time-stepping
        particle.property<position>() += globalVar.property<dt>() * particle.property<velocity>();
    }

    void interact(GlobalVar<globalvar_type> globalVar, Particle<dimension, position_type, property_type_n<dimension>> particle,
                  Particle<dimension, position_type, property_type_n<dimension>> neighbor) override {

        // Check cutoff radius
        auto diff = particle.property<position>() - neighbor.property<position>();
        if (diff > globalVar.property<r_cut>())
            return;

        // Compute projection of collision
        auto diff_norm = diff / diff.abs2();
        auto collision_p = (diff * particle.property<velocity>()) * diff_norm;
        auto collision_n = (diff * neighbor.property<velocity>()) * diff_norm;
        auto diff_collision = collision_n - collision_p;

        // Apply collision
        particle.property<velocity>() += diff_collision;
        neighbor.property<velocity>() -= diff_collision;
    }

    void evolveGlobalVar(GlobalVar<globalvar_type> globalVar) override {
        globalVar.property<time>() += globalVar.property<dt>();
    }

    void stop(GlobalVar<globalvar_type> globalVar) override {
        return globalVar.property<time>() >= globalVar.property<t_final>();
    }

};

#endif //OPENFPM_PDATA_DEM_TEST_HPP
