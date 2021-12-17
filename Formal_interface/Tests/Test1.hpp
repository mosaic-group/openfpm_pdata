//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_TEST1_HPP
#define OPENFPM_PDATA_TEST1_HPP


#include <Vector/vector_dist.hpp>
#include "../Particle.hpp"
#include "../ParticleData.hpp"
#include "../ParticleMethod.hpp"
#include "../Transition.hpp"

//typedef aggregate<float, float> particle_type;
//typedef aggregate<float, float, float, float> globalvar_type;
//Property<velocity>(particle)

// Position type
using position_type = float;

// Particle type
template <int dimension>
using particle_type_n = aggregate<float[dimension], float[dimension]>;


template <int dimension>
class Test1 : public ParticleMethod<dimension, position_type, particle_type_n<dimension>> {
public:

    bool freeParticles;

    constexpr static position_type domainMin = 0.0;
    constexpr static position_type domainMax = 1.0;

    typedef Particle<typename Test1::particleType> Particle_;

    static constexpr int position = 0;
    static constexpr int velocity = 1;

    static constexpr int dt = 0;
    static constexpr int time = 1;
    static constexpr int t_final = 2;
    static constexpr int r_cut = 3;

    int iteration = 0;

    void evolve(ParticleRef<dimension, position_type, particle_type_n<dimension>> particle) override {

//        std::cout << "evolution " << particle.template property<position>()[0] << " , " << particle.template property<velocity>()[0] << std::endl;

        // Euler time-stepping
        particle.template property<position>()[0] += particle.template property<velocity>()[0];

//        std::cout << "evolution " << particle.template property<position>()[0] << " , " << particle.template property<velocity>()[0] << std::endl;

    }
/*
    void interact(Particle<particle_type_n<dimension>> particle, Particle<particle_type_n<dimension>> neighbor) override {

    }*/

    bool stop() override {
        iteration++;

//        std::cout << "stop it " << iteration << std::endl;
        if (iteration > 2)
            return true;
        return false;
    }
};


#endif //OPENFPM_PDATA_TEST1_HPP
