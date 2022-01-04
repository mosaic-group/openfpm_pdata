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
#include "../InitialCondition.hpp"

//typedef aggregate<float, float> particle_type;
//typedef aggregate<float, float, float, float> globalvar_type;
//Property<velocity>(particle)

// Position type
using position_type = float;

// Particle type
template <int dimension>
using particle_type_n = aggregate<float[dimension], float[dimension]>;

// GlobalVar type
using globalvar_type = aggregate<float, float, float, float>;

template <int dimension>
class Test1 : public ParticleMethod<dimension, position_type, particle_type_n<dimension>, globalvar_type> {
public:

    struct GV {
        float dt = 0.1;
        float t = 0;
        float t_final = 30;
        float r_cut = 0.5;
    } globalvar;

    bool freeParticles;

    constexpr static position_type domainMin = 0.0;
    constexpr static position_type domainMax = 10.0;

//    typedef Particle<typename Test1::particleType> Particle_;

    static constexpr int position = 0;
    static constexpr int velocity = 1;

    static constexpr int dt = 0;
    static constexpr int time = 1;
    static constexpr int t_final = 2;
    static constexpr int r_cut = 3;

    int iteration = 0;

    void evolve(/*GlobalVar<globalvar_type> globalVar,*/ Particle<dimension, position_type, particle_type_n<dimension>> particle) override {

//        std::cout << "evolution " << particle.template property<position>()[0] << " , " << particle.template property<velocity>()[0] << std::endl;

        // Euler time-stepping
        particle.template position()[0] += particle.template property<velocity>()[0] * globalvar.dt;
        particle.template position()[1] += particle.template property<velocity>()[1] * globalvar.dt;

        //        particle.template property<position>()[0] += particle.template property<velocity>()[0] * globalvar.dt;

//        std::cout << "evolution " << particle.template property<position>()[0] << " , " << particle.template property<velocity>()[0] << std::endl;

    }
/*
    void interact(Particle<particle_type_n<dimension>> particle, Particle<particle_type_n<dimension>> neighbor) override {

    }*/

    bool stop() override {
/*        iteration++;

//        std::cout << "stop it " << iteration << std::endl;
        if (iteration > 2)
            return true;
        return false;*/

        if (globalvar.t > globalvar.t_final)
            return true;
        globalvar.t += globalvar.dt;
        return false;
    }
};

template <typename ParticleMethodType>
class InitialCondition1 : InitialCondition<ParticleMethodType> {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

public:
    constexpr static PositionType domainMin[dimension] = {0.0, 0.0};
    constexpr static PositionType domainMax[dimension] = {10.0, 10.0};
    constexpr static size_t boundaryCondition = PERIODIC;

};



#endif //OPENFPM_PDATA_TEST1_HPP
