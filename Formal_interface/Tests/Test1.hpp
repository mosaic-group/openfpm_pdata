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
#include "../Neighborhood.hpp"


struct PSE_ParticleSignature {
    static constexpr int dimension = 2;
    typedef double position;
    typedef aggregate<float, float> properties;
    typedef MESH_PARTICLES dataStructure;
};

// Property identifier
constexpr int concentration = 0;
constexpr int accumulator = 1;




struct GlobalVariable {
    float dt = 0.05;
    float t = 0;
    float t_final = 10.1;

    float domainSize = 40.0;
    int meshSize = 256;
    float meshSpacing = domainSize / meshSize;
    float epsilon = meshSpacing;
    float r_cut = 3 * epsilon;
    float D = 0.01;
    float kernel = dt * D * 15.0 * pow(meshSpacing/epsilon, 3)  / pow(epsilon * M_PI, 2);
} globalvar;




// Particle Method implementation

template <typename ParticleSignature>
class PSE_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {
        Point<dimension, PositionType> p_pos = particle.position();
        Point<dimension, PositionType> n_pos = neighbor.position();
        PositionType distance2 = p_pos.distance2(n_pos);

        particle.template property<accumulator>() +=
                (neighbor.template property<concentration>() - particle.template property<concentration>())
                / (1 + pow(distance2 / globalvar.epsilon / globalvar.epsilon, 5)) ;

    }

    void evolve(Particle<ParticleSignature> particle) override {

            particle.template property<concentration>() += globalvar.kernel * particle.template property<accumulator>();
            particle.template property<accumulator>()= 0;

    }

    bool stop() override {
        std::cout << "\r" << int(globalvar.t / globalvar.t_final * 100) << "%" << std::flush;

        // Check simulation time
        if (globalvar.t > globalvar.t_final) {
            std::cout << std::endl;
            return true;
        }

        globalvar.t += globalvar.dt;

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

    void initialization(Particle<ParticleSignatureType> particle) override {

        // Randomize concentration (normal distribution)
        particle.template property<concentration>() = this->normalDistribution(0, 5);
    }

};



#endif //OPENFPM_PDATA_TEST1_HPP
