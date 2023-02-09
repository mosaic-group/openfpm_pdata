//
// Created by landfried on 27.09.22.
//

#ifndef OPENFPM_PDATA_EC_ALGORITHM_HPP
#define OPENFPM_PDATA_EC_ALGORITHM_HPP


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
#include <random>
#include <time.h>


struct DEM_ParticleSignature {
    static constexpr int dimension = 2;
    typedef float position;
    typedef aggregate<float[dimension], float[dimension]> properties;
    typedef FREE_PARTICLES dataStructure;
};

// Property identifier
constexpr int velocity = 0;
constexpr int deltaVelocity = 1;

struct GlobalVariable {
    static float dt;
    static float t;
    static float t_final;
    static float r_cut;
    static float domainSize;
} globalvar;

class DEM_ParticleMethod : public ParticleMethod<DEM_ParticleSignature> {

public:
    // std::random_device rd;
    
    // static std::random_device dev;
    // static std::mt19937 rng();
    void evolve(Particle<ParticleSignature> particle) override {
        // std::mt19937 rng;
        // std::normal_distribution<> dis(0.0,M_PI*0.05);
        // srand(time(NULL));
        float rx=((float)rand() / (float)RAND_MAX - .5)/8.0;
        float ry=((float)rand() / (float)RAND_MAX - .5)/8.0;
        float rz=((float)rand() / (float)RAND_MAX - .5)/8.0;
        // float rx=dis(rng);
        // dis.reset();
        // float ry=dis(rng);
        // dis.reset();
        // float rz=dis(rng);
        PARTICLE(deltaVelocity) += PARTICLE(velocity)*1000.0; 
        Point<dimension,float> a = particle.template property_raw<deltaVelocity>();
        // std::cout << a.toString() << std::endl;
        // std::cout << PARTICLE(deltaVelocity)[0]*PARTICLE(deltaVelocity)[0] <<", "<< PARTICLE(deltaVelocity)[1]*PARTICLE(deltaVelocity)[1] <<", "<< PARTICLE(deltaVelocity)[2]*PARTICLE(deltaVelocity)[2] << std::endl;
        float absolute = a[0]*a[0]
                        +a[1]*a[1];
                        // +a[2]*a[2];
        // std::cout << absolute << std::endl;
        // // Apply change of velocity
        
        PARTICLE(velocity) = PARTICLE(deltaVelocity)/sqrt(absolute);
        // // std::cout << PARTICLE(velocity)[0] <<", "<< PARTICLE(velocity)[1] <<", "<< PARTICLE(velocity)[2] << std::endl;

        // PARTICLE(velocity)[0]=   cos(ry)*cos(rz)*PARTICLE(velocity)[0]
        //                         +(sin(rx)*sin(ry)*cos(rz)-cos(rx)*sin(rz))*PARTICLE(velocity)[1]
        //                         +(cos(rx)*sin(ry)*cos(rz)+sin(rx)*sin(rz))*PARTICLE(velocity)[2];
        
        // PARTICLE(velocity)[1]=   cos(ry)*sin(rz)*PARTICLE(velocity)[0]
        //                         +(sin(rx)*sin(ry)*sin(rz)+cos(rx)*cos(rz))*PARTICLE(velocity)[1]
        //                         +(cos(rx)*sin(ry)*sin(rz)-sin(rx)*cos(rz))*PARTICLE(velocity)[2];

        // PARTICLE(velocity)[2]=  -sin(ry)*PARTICLE(velocity)[0]
        //                         +sin(rx)*cos(ry)*PARTICLE(velocity)[1]
        //                         +cos(rx)*cos(ry)*PARTICLE(velocity)[2];                      

        PARTICLE(velocity)[0] += rx;
        PARTICLE(velocity)[1] += ry;
        // PARTICLE(velocity)[2] += rz;
        
        // normailze
        Point<dimension,float> magnitude = particle.template property_raw<velocity>();
        float absolute2 = magnitude[0]*magnitude[0]
                        +magnitude[1]*magnitude[1];
                        // +magnitude[2]*magnitude[2];

        PARTICLE(velocity) /= sqrt(absolute2);
        
        //Euler time-stepping move particles
        particle.position() += PARTICLE(velocity) * globalvar.dt;

        //std::cout << PARTICLE(deltaVelocity)[0] << std::endl;
        // Reset change of velocity
        PARTICLE(deltaVelocity) = 0.0f;
    }

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {
        // Compute collision
        PARTICLE(deltaVelocity) +=NEIGHBOR(velocity);
        //PARTICLE(deltaVelocity)[0]++;
        // std::cout << "interact" << std::endl;
        // std::cout << PARTICLE(velocity)[0] <<", "<< PARTICLE(velocity)[1] <<", "<< PARTICLE(velocity)[2] << std::endl;
    }


    void evolveGlobalVariable() {
        // advance time
        globalvar.t += globalvar.dt;
    }

    bool stop() override {
        // Check simulation time
        return globalvar.t > globalvar.t_final;
    }

};



#endif //OPENFPM_PDATA_EC_ALGORITHM_HPP
