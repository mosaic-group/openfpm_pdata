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
#include "../SimulationParameters.hpp"
#include "../InitialCondition.hpp"
#include "../Neighborhood.hpp"




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
class DEM_ParticleMethod : public ParticleMethod<dimension, position_type, property_type_n<dimension>, globalvar_type> {

public:

    struct GlobalVariable {
        float dt = 0.04;
        float t = 0;
        float t_final = 50;
        float r_cut = 0.3;
        float damp = 0.9;
    } globalvar;

    void evolve(Particle<dimension, position_type, property_type_n<dimension>> particle) override {

        for (int i = 0; i < dimension; i++) {

            // Apply change of velocity
            particle.template property<velocity>()[i] += particle.template property<acceleration>()[i];

            // Reset change of velocity
            particle.template property<acceleration>()[i] = 0.0f;

            // Euler time-stepping move particles
            particle.template position()[i] += particle.template property<velocity>()[i] * globalvar.dt;
        }

    }

    void interact(Particle<dimension, position_type, property_type_n<dimension>> particle, Particle<dimension, position_type, property_type_n<dimension>> neighbor) override {

        // Declare particle property variables
        Point<dimension, position_type> p_pos = particle.position();
        Point<dimension, position_type> n_pos = neighbor.position();
        Point<dimension, position_type> p_vel = particle.template property<velocity>();
        Point<dimension, position_type> n_vel = neighbor.template property<velocity>();


        // Check cutoff radius
        if (p_pos.distance(n_pos) > globalvar.r_cut)
            return;

        // Check if particles are moving towards each other
        Point<dimension, position_type> p_move = p_pos + (p_vel * globalvar.dt);
        Point<dimension, position_type> n_move = n_pos + (n_vel * globalvar.dt);
        float dist = p_pos.distance2(n_pos);
        float dist_move = (p_move).distance2(n_move);
        if (dist < dist_move)
            return;

        // Compute collision vector
        Point<dimension, position_type> diff = p_pos - n_pos;
        Point<dimension, position_type> diff_scaled = diff / p_pos.distance2(n_pos);
        Point<dimension, position_type> p_collision = (diff * p_vel) * diff_scaled;
        Point<dimension, position_type> n_collision = (diff * n_vel) * diff_scaled;
        Point<dimension, position_type> diff_collision = n_collision - p_collision;

        diff_collision = diff_collision * globalvar.damp;

        // Apply collision to particle acceleration

        for (int i = 0; i < dimension; i++) {
            particle.template property<acceleration>()[i] += diff_collision[i];
        }
    }


    bool stop() override {
        std::cout << globalvar.t << std::endl;

        // Check simulation time
        if (globalvar.t > globalvar.t_final)
            return true;

        globalvar.t += globalvar.dt;

        return false;
    }
};

template <typename ParticleMethodType>
class DEM_SimulationParams : public SimulationParameters<ParticleMethodType> {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    ParticleMethodType particleMethod;

public:

    // Domain
//    PositionType domainMin[dimension] = {0.0, 0.0};
//    PositionType domainMax[dimension] = {20.0, 20.0};

    Point<dimension, PositionType> domainMin;
    Point<dimension, PositionType> domainMax;

//    size_t meshSize[dimension] = {10, 10, 10};


    DEM_SimulationParams() : domainMin(0.0f), domainMax(20.0f) {}


//    Point<dimension, PositionType> (0.0);

    // Boundary conditions
    size_t boundaryConditions[dimension] = {PERIODIC, PERIODIC};

/*
    // Mesh initial condition
    typedef InitialConditionMesh initialCondition;
    constexpr static size_t meshSize[dimension] = {18, 18};
*/

    // Random initial condition
    typedef InitialConditionRandom initialCondition;
    int numberParticles = 300;

    // Neighborhood method
    typedef CellListNeighborhood neighborhoodDetermination;
    float cellWidth = particleMethod.globalvar.r_cut;

    void initialization(Particle<dimension, PositionType , PropertyType> particle) override {

        // Randomize velocity (normal distribution)
        for (int i = 0; i < dimension; i++) {
            particle.template property<velocity>()[i] = this->normalDistribution(0, 2);
        }
    }

};


#endif //OPENFPM_PDATA_DEM_TEST_HPP
