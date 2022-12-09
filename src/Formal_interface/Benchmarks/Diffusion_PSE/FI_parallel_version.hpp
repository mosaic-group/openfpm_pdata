//
// Created by landfried on 17.12.21.
//

#ifndef OPENFPM_PDATA_PSE_TEST_HPP
#define OPENFPM_PDATA_PSE_TEST_HPP

#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"
#include "Formal_interface/Instance.hpp"


#define PARTICLE(property_arg) particle.template property<property_arg>()
#define NEIGHBOR(property_arg) neighbor.template property<property_arg>()

struct PSE_ParticleSignature {
    static constexpr int dimension = 3;
    typedef double position;
    typedef aggregate<size_t, double, double> properties;
    typedef MESH_PARTICLES dataStructure;
};

// Property identifier
constexpr int concentration = 1;
constexpr int accumulator = 2;




struct GlobalVariable {
    double dt = 0.05;
    double t = 0;
    double t_final = 0.5;

    double domainSize = 40.0;
    int meshSize = 128;
    double meshSpacing = domainSize / meshSize;
    double epsilon = meshSpacing;
    double r_cut = 3 * epsilon;
    double D = 0.01;
    double kernel = dt * D * 15.0 * pow(meshSpacing/epsilon, PSE_ParticleSignature::dimension  )  / pow(epsilon * M_PI, 2);
} globalvar;




// Particle Method implementation

template <typename ParticleSignature>
class PSE_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {
        Point<dimension, PositionType> p_pos = particle.position_raw();
        Point<dimension, PositionType> n_pos = neighbor.position_raw();
        PositionType distance2 = p_pos.distance2(n_pos);

        double exchange = (NEIGHBOR(concentration) - PARTICLE(concentration))
                                                  / (1 + pow(distance2 / globalvar.epsilon / globalvar.epsilon, 5)) ;
// cottet kernel

//        particle.template property<accumulator>() += exchange;
//        neighbor.template property<accumulator>() -= exchange;
        PARTICLE(accumulator) += exchange;
        NEIGHBOR(accumulator) -= exchange;

    }

    void evolve(Particle<ParticleSignature> particle) override {

//        particle.template property<concentration>() += particle.template property<accumulator>() * globalvar.kernel;
//        particle.template property<accumulator>() = 0;
        PARTICLE(concentration) += PARTICLE(accumulator) * globalvar.kernel;
        PARTICLE(accumulator) = 0;

    }

    void evolveGlobalVariable() override {

        // advance time
        globalvar.t += globalvar.dt;
//        std::cout << "\r" << int(globalvar.t / globalvar.t_final * 100) << "%" << std::flush;

    }

    bool stop() override {

        // Check simulation time
        if (globalvar.t > globalvar.t_final)
            return true;

        return false;
    }
};

template <typename ParticleSignatureType>
class PSE_SimulationParams : public SimulationParameters<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    PSE_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
        this->setMeshSize(globalvar.meshSize);
        this->setCutoffRadius(globalvar.r_cut);
        this->setMeshSpacing(globalvar.meshSpacing);
    }

    typedef NEIGHBORHOOD_MESH neighborhoodDetermination;
    static const int interactionType = INTERACTION_SYMMETRIC;

//    int writeIteration = 100;
        bool writeOutput = false;


    void initialization(Particle<ParticleSignatureType> particle) override {

/*        // Randomize concentration (normal distribution)
        particle.template property<concentration>() = this->normalDistribution(0, 5);*/

        bool centerParticle = true;
        for (int i = 0; i < dimension; ++i) {
            if (particle.position()[i] != globalvar.domainSize / 2)
                centerParticle = false;
        }

        if (centerParticle) {
            PARTICLE(concentration) = 1;
        }

    }

};



#endif //OPENFPM_PDATA_PSE_TEST_HPP
