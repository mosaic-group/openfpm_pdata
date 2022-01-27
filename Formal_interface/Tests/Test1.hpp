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
    typedef aggregate<double, double> properties;
    typedef MESH_PARTICLES dataStructure;
};

// Property identifier
constexpr int concentration = 0;
constexpr int accumulator = 1;

template <typename ParticleSignature>
class PSE_ParticleMethod : public ParticleMethod<ParticleSignature> {

    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;

public:

    struct GlobalVariable {
        float dt = 0.04;
        float t = 0;
        float t_final = 20.1;

        float domainSize = 2.5;
        int meshSize = 128;
        float meshSpacing = domainSize / meshSize;
        float epsilon = meshSpacing;
        float r_cut = 3 * epsilon;
        float D = 0.01;
        float kernel = dt * D * 15.0 * ((meshSpacing / epsilon) * (meshSpacing / epsilon) * (meshSpacing / epsilon)) / ((epsilon * M_PI) * (epsilon * M_PI));
    } globalvar;




    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {

        Point<dimension, PositionType> p_pos = particle.position();
        Point<dimension, PositionType> n_pos = neighbor.position();
        PositionType distance2 = p_pos.distance2(n_pos);

        particle.template property<accumulator>() +=
                (neighbor.template property<concentration>() - particle.template property<concentration>())
                / (((distance2 / globalvar.epsilon / globalvar.epsilon) *
                    (distance2 / globalvar.epsilon / globalvar.epsilon) *
                    (distance2 / globalvar.epsilon / globalvar.epsilon) *
                    (distance2 / globalvar.epsilon / globalvar.epsilon) *
                    (distance2 / globalvar.epsilon / globalvar.epsilon))
                   + 1.0);

    }

    void evolve(Particle<ParticleSignature> particle) override {

            particle.template property<concentration>() += globalvar.kernel * particle.template property<accumulator>();
            particle.template property<accumulator>()= 0;

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

template <typename ParticleSignatureType>
class PSE_SimulationParams : public SimulationParameters<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    // Domain
    Point<dimension, PositionType> domainMin;
    Point<dimension, PositionType> domainMax;

    // Boundary conditions
    periodicity<dimension> boundaryConditions {PERIODIC, PERIODIC};

    size_t meshSize[dimension] = {128, 128};

    PSE_SimulationParams() : domainMin(0.0f), domainMax(2.5f) {
//        std::fill(std::begin(boundaryConditions), std::end(boundaryConditions), PERIODIC);
    }

/*
    // Mesh initial condition
    typedef InitialConditionMesh initialCondition;
    constexpr static size_t meshSize[dimension] = {18, 18};
*/

    typedef MeshNeighborhood neighborhoodDetermination;

    void initialization(Particle<ParticleSignatureType> particle) override {

        // Randomize velocity (normal distribution)
        particle.template property<concentration>() = this->normalDistribution(0, 2);
//            particle.template property<concentration>()[i] = particle.position()[i];
    }

};



#endif //OPENFPM_PDATA_TEST1_HPP
