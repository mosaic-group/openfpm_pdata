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
#include "../SimulationParameters.hpp"
#include "../InitialCondition.hpp"
#include <valarray>

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
        float dt = 0.05;
        float t = 0;
        float t_final = 0.1;
        float r_cut = 0.3;
    } globalvar;

    constexpr static position_type domainMin = 0.0;
    constexpr static position_type domainMax = 20.0;

    static constexpr int velocity = 0;
    static constexpr int acceleration = 1;

    void evolve(/*GlobalVar<globalvar_type> globalVar,*/ Particle<dimension, position_type, particle_type_n<dimension>> particle) override {

//        std::cout << "evolution " << particle.template property<position>()[0] << " , " << particle.template property<velocity>()[0] << std::endl;

        // Apply change of velocity
        particle.template property<velocity>()[0] += particle.template property<acceleration>()[0];
        particle.template property<velocity>()[1] += particle.template property<acceleration>()[1];

        // Reset change of velocity
        particle.template property<acceleration>()[0] = 0.0f;
        particle.template property<acceleration>()[1] = 0.0f;

        // Euler time-stepping move particles
        particle.template position()[0] += particle.template property<velocity>()[0] * globalvar.dt;
        particle.template position()[1] += particle.template property<velocity>()[1] * globalvar.dt;

        //        particle.template property<position>()[0] += particle.template property<velocity>()[0] * globalvar.dt;

//        std::cout << "evolution " << particle.template property<position>()[0] << " , " << particle.template property<velocity>()[0] << std::endl;

    }

    void interact(Particle<dimension, position_type, particle_type_n<dimension>> particle, Particle<dimension, position_type, particle_type_n<dimension>> neighbor) override {
//        std::cout << "interact" << std::endl;

        Point<dimension, position_type> p_pos = particle.position();
        Point<dimension, position_type> n_pos = neighbor.position();
        Point<dimension, position_type> p_vel = particle.template property<velocity>();
        Point<dimension, position_type> n_vel = neighbor.template property<velocity>();

        if (p_pos.distance(n_pos) > globalvar.r_cut)
            return;

        Point<dimension, position_type> diff = p_pos - n_pos;
        Point<dimension, position_type> diff_scaled = diff / p_pos.distance2(n_pos);
        Point<dimension, position_type> p_collision = (diff * p_vel) * diff_scaled;
        Point<dimension, position_type> n_collision = (diff * n_vel) * diff_scaled;
        Point<dimension, position_type> diff_collision = n_collision - p_collision;

//        p_vel += diff_collision;
//        n_vel -= diff_collision;

        particle.template property<acceleration>()[0] = diff_collision.asArray()[0];
        particle.template property<acceleration>()[1] = diff_collision.asArray()[1];

//        particle.template property<velocity>()[0] = p_vel.asArray()[0];
//        particle.template property<velocity>()[1] = p_vel.asArray()[1];
//        neighbor.template property<velocity>() = n_vel.asArray();

    }


    bool stop() override {
/*        iteration++;

//        std::cout << "stop it " << iteration << std::endl;
        if (iteration > 2)
            return true;
        return false;*/

        std::cout << globalvar.t << std::endl;

        if (globalvar.t > globalvar.t_final)
            return true;
        globalvar.t += globalvar.dt;
        return false;
    }
};

template <typename ParticleMethodType>
class SimulationParams1 : public SimulationParameters<ParticleMethodType> {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

public:

    // Domain
//    constexpr static PositionType domainMin[dimension] = {0.0, 0.0};
//    constexpr static PositionType domainMax[dimension] = {20.0, 20.0};
    Point<dimension, PositionType> domainMin = {0.0, 0.0};
    Point<dimension, PositionType> domainMax = {20.0, 20.0};
//    PositionType domainMin[dimension] = {0.0, 0.0};
//    PositionType domainMax[dimension] = {20.0, 20.0};

    // Boundary conditions
    size_t boundaryCondition = PERIODIC;

/*
    // Mesh initial condition
    typedef InitialConditionMesh initialCondition;
    constexpr static size_t meshSize[dimension] = {5, 5};
*/

    // Random initial condition
    typedef InitialConditionRandom initialCondition;
    int numberParticles = 50;


};



#endif //OPENFPM_PDATA_TEST1_HPP
