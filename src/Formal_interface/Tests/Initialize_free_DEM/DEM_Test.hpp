//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_DEM_TEST_HPP
#define OPENFPM_PDATA_DEM_TEST_HPP


#include <array>
#include <Vector/vector_dist.hpp>
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"
#include "Formal_interface/Instance.hpp"


struct DEM_ParticleSignature {
    static constexpr int dimension = 2;
    typedef float position;
    typedef aggregate<float[dimension], float[dimension]> properties;
    typedef FREE_PARTICLES dataStructure;
};

// Property identifier
constexpr int velocity = 0;
constexpr int acceleration = 1;


struct GlobalVariable {
    float dt = 0.04;
    float t = 0;
    float t_final = 20.1;
    float r_cut = 0.5;
    float damp = 0.9;
    float domainSize = 20.0;
} globalvar;


template <typename ParticleSignature>
class DEM_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:


    void evolve(Particle<ParticleSignature> particle) override {

        // Apply change of velocity
        particle.template property_vec<velocity>() += particle.template property_vec<acceleration>();

        // Reset change of velocity
        particle.template property_vec<acceleration>() = 0.0f;

        // Euler time-stepping move particles
        particle.position_vec() += particle.template property_vec<velocity>() * globalvar.dt;

    }

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {

//        std::cout << "interact" << std::endl;

        // Declare particle property variables
        Point<dimension, PositionType> p_pos = particle.position();
        Point<dimension, PositionType> n_pos = neighbor.position();
        Point<dimension, PositionType> p_vel = particle.template property<velocity>();
        Point<dimension, PositionType> n_vel = neighbor.template property<velocity>();


        // Check cutoff radius
        if (p_pos.distance(n_pos) > globalvar.r_cut)
            return;

        // Check if particles are moving towards each other
        Point<dimension, PositionType> p_move = p_pos + (p_vel * globalvar.dt);
        Point<dimension, PositionType> n_move = n_pos + (n_vel * globalvar.dt);
        float dist = p_pos.distance2(n_pos);
        float dist_move = (p_move).distance2(n_move);
        if (dist < dist_move)
            return;

        // Compute collision vector
        Point<dimension, PositionType> diff = p_pos - n_pos;
        Point<dimension, PositionType> diff_scaled = diff / p_pos.distance2(n_pos);
        Point<dimension, PositionType> p_collision = (diff * p_vel) * diff_scaled;
        Point<dimension, PositionType> n_collision = (diff * n_vel) * diff_scaled;
        Point<dimension, PositionType> diff_collision = n_collision - p_collision;

        diff_collision = diff_collision * globalvar.damp;

        // Apply collision to particle acceleration
        particle.template property_vec<acceleration>() += diff_collision;
//        neighbor.template property_vec<acceleration>() -= diff_collision;

    }

    void evolveGlobalVariable() override {

        // advance time
        globalvar.t += globalvar.dt;
        std::cout << "\r" << int(globalvar.t / globalvar.t_final * 100) << "%" << std::flush;

    }

    bool stop() override {

        // Check simulation time
        if (globalvar.t > globalvar.t_final)
            return true;

        return false;
    }
};

template <typename ParticleSignatureType>
class DEM_SimulationParams : public SimulationParameters<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:


    DEM_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
        this->setCutoffRadius(globalvar.r_cut);
        this->setNumberParticles(0);
        this->setCellWidth(globalvar.r_cut);
    }

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;

    bool writeOutput = true;

/*    void initialization(Particle<ParticleSignatureType> particle) override {

        // Randomize velocity (normal distribution)
        for (int i = 0; i < dimension; i++) {
            particle.template property<velocity>()[i] = this->normalDistribution(0, 3);
        }
    }*/

};

class Instance1 : Instance<DEM_ParticleMethod<DEM_ParticleSignature>, DEM_SimulationParams<DEM_ParticleSignature>> {

public:

    Instance1(ParticleData<DEM_ParticleMethod<DEM_ParticleSignature>, DEM_SimulationParams<DEM_ParticleSignature>> &particleData_in) :
        Instance<DEM_ParticleMethod<DEM_ParticleSignature>, DEM_SimulationParams<DEM_ParticleSignature>>(particleData_in){}

    virtual void freePlacement() {

        // block
        for (int i = 4; i < 18; ++i) {
            for (int j = 1; j < 3; ++j) {
                this->addParticle();
                this->position()[0] = float(i);
                this->position()[1] = float(j);
                this->property<velocity>()[0] = 0.0;
                this->property<velocity>()[1] = 0.0;
            }
        }

        // cirlce
        for (int i = 0; i < 30; ++i) {
            this->addParticle();
            this->position()[0] = sin(M_PI * float(i) / 15) * 5 + 10;
            this->position()[1] = cos(M_PI * float(i) / 15) * 5 + 12;
            this->property<velocity>()[0] = -sin(M_PI * float(i) / 15) * float(i) / 7;
            this->property<velocity>()[1] = -cos(M_PI * float(i) / 15) * float(i) / 7;
        }

    }
};


#endif //OPENFPM_PDATA_DEM_TEST_HPP
