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



// Position type
using position_type = float;

// Property type
template <int dimension>
using property_type_n = aggregate<float[dimension], float[dimension]>;

// Property identifier
constexpr int velocity = 0;
constexpr int acceleration = 1;

// GlobalVar type
using globalvar_type = aggregate<float, float, float, float>;

template <int dimension>
class Test1 : public ParticleMethod<dimension, position_type, property_type_n<dimension>, globalvar_type> {

public:

    struct GlobalVariable {
        float dt = 0.05;
        float t = 0;
        float t_final = 10;
        float r_cut = 0.3;
    } globalvar;




    void initialization(Particle<dimension, position_type, property_type_n<dimension>> particle) override {

        for (int i = 0; i < dimension; i++) {
            particle.template property<velocity>()[i] = this->normalDistribution(0, .5);
        }

    }

    void evolve(/*GlobalVar<globalvar_type> globalVar,*/ Particle<dimension, position_type, property_type_n<dimension>> particle) override {

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

    void interact(Particle<dimension, position_type, property_type_n<dimension>> particle, Particle<dimension, position_type, property_type_n<dimension>> neighbor) override {
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

        particle.template property<acceleration>()[0] = diff_collision[0];
        particle.template property<acceleration>()[1] = diff_collision[1];
    }


    bool stop() override {
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
    PositionType domainMin[dimension] = {0.0, 0.0};
    PositionType domainMax[dimension] = {20.0, 20.0};

    // Boundary conditions
    size_t boundaryConditions[dimension] = {PERIODIC, PERIODIC};

/*
    // Mesh initial condition
    typedef InitialConditionMesh initialCondition;
    constexpr static size_t meshSize[dimension] = {5, 5};
*/

    // Random initial condition
    typedef InitialConditionRandom initialCondition;
    int numberParticles = 50;



/*
    void initialConditions() {
        create_random(50);
        create_mesh({5, 5});
        random(velocity, dimension, -.4, .4);
        random(acceleration, dimension, -.4, .4);
    }

    std::map<int, InitialCondition> PropertyInitialCondition;
    PropertyInitialCondition[velocity] =

    InitialConditionPosition initialConditionPosition(random_uniform, 0.0, 20.0);
    InitialConditionProperty initialConditionVelocity()
    InitialCondition(position, random_uniform, 0.0, 20.0);
    InitialCondition(position, mesh, 0.0, 20.0);
    InitialCondition(velocity, random, 0.0, 20.0);
*/

};



#endif //OPENFPM_PDATA_TEST1_HPP
