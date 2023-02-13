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
    static bool addparticles;
} globalvar;



class DEM_ParticleMethod : public ParticleMethod<DEM_ParticleSignature> {

public:

    void evolve(Particle<ParticleSignature> particle) override {
        // Apply change of velocity
        PARTICLE_WRITE(velocity) += PARTICLE(deltaVelocity);
        // Reset change of velocity
        PARTICLE_WRITE(deltaVelocity) = 0.0f;
        // Euler time-stepping move particles
        PARTICLE_POSITION_WRITE += PARTICLE(velocity) * globalvar.dt;
    }

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {
        // Compute collision
        auto diff = NEIGHBOR_POSITION - PARTICLE_POSITION;
        PARTICLE_WRITE(deltaVelocity) += diff / norm2(diff) * pmul(diff, NEIGHBOR(velocity) - PARTICLE(velocity));
    }

   // cmake ../.  -DBOOST_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/BOOST -DHDF5_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/HDF5/ -DLIBHILBERT_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/LIBHILBERT -DMPI_VENDOR=openmpi -DPARMETIS_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/PARMETIS -DMETIS_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/METIS -DPETSC_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/PETSC -DSUITESPARSE_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/SUITESPARSE -DEIGEN3_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/EIGEN -DOPENBLAS_ROOT=/home/landfried/Studium/OpenFPM/2021_9_17/dependencies/OPENBLAS/ -DCMAKE_BUILD_TYPE=Release



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
