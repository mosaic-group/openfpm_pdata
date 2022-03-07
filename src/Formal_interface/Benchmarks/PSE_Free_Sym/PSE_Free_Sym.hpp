//
// Created by landfried on 02.03.22.
//

#ifndef OPENFPM_PDATA_PSE_FREE_UNSYM_HPP
#define OPENFPM_PDATA_PSE_FREE_UNSYM_HPP


#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Neighborhood.hpp"


struct PSE_ParticleSignature {
    static constexpr int dimension = 2;
    typedef float position;
    typedef aggregate<float, float> properties;
    typedef FREE_PARTICLES dataStructure;
};

// Property identifier
constexpr int concentration = 0;
constexpr int accumulator = 1;

struct GlobalVariable {
    float dt = 0.05;
    float t = 0;
    float t_final = 10.1;

    float domainSize = 200.0;
    size_t meshSize = 256;
    float meshSpacing = domainSize / (float)meshSize;
    float epsilon = meshSpacing;
    float r_cut = 3 * epsilon;
    float D = 0.01;
    float kernel = dt * D * 15.0 * pow(meshSpacing/epsilon, 3)  / pow(epsilon * M_PI, 2);

    void recalculate() {
        meshSpacing = domainSize / (float)meshSize;
        epsilon = meshSpacing;
        r_cut = 3 * epsilon;
        kernel = dt * D * 15.0 * pow(meshSpacing/epsilon, 3)  / pow(epsilon * M_PI, 2);
    }

} globalvar;


template <typename ParticleSignature>
class Benchmark_ParticleMethod1 : public ParticleMethod<ParticleSignature> {
    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:


    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {
        Point<dimension, PositionType> p_pos = particle.position();
        Point<dimension, PositionType> n_pos = neighbor.position();

//        std::cout << p_pos.toString() << " <-> " << n_pos.toString() << std::endl;

        PositionType distance2 = p_pos.distance2(n_pos);

        float exchange = (neighbor.template property_vec<concentration>() - particle.template property_vec<concentration>())
                         / (1 + pow(distance2 / globalvar.epsilon / globalvar.epsilon, 5)) ;

        particle.template property_vec<accumulator>() += exchange;
        neighbor.template property_vec<accumulator>() -= exchange;

    }

    void evolve(Particle<ParticleSignature> particle) override {

        particle.template property_vec<concentration>() += particle.template property_vec<accumulator>() * globalvar.kernel;
        particle.template property_vec<accumulator>() = 0;

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
class Benchmark_SimulationParams : public SimulationParameters<ParticleSignatureType> {
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    Benchmark_SimulationParams() {
        this->setDomain(globalvar.domainSize);
        this->setBoundaryConditions(PERIODIC);
        this->setCutoffRadius(globalvar.r_cut);
        this->setMeshSize(globalvar.meshSize);
        this->setCellWidth(globalvar.r_cut);
    }

    // Mesh initial condition
    typedef INITIALCONDITION_MESH initialCondition;

    // Neighborhood method
    typedef NEIGHBORHHOD_CELLLIST neighborhoodDetermination;

    static const int interactionType = INTERACTION_SYMMETRIC;

    // Output
    bool writeOutput = false;

    void initialization(Particle<ParticleSignatureType> particle) override {

//        std::cout << "r_cut " << globalvar.r_cut << std::endl;
        for (int i = 0; i < dimension; i++) {
            // Randomize concentration (normal distribution)
            particle.template property<concentration>() = this->normalDistribution(0, 5);
        }
    }



};



#endif //OPENFPM_PDATA_PSE_FREE_UNSYM_HPP
