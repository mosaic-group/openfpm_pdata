//
// Created by landfried on 29.09.22.
//

#ifndef OPENFPM_PDATA_TP_ALGORITHM_HPP
#define OPENFPM_PDATA_TP_ALGORITHM_HPP

#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Interaction_Impl.hpp"




struct PSE_ParticleSignature {
    static constexpr int dimension = 2;
    typedef double position;
    typedef aggregate<size_t, double, double> properties;
    typedef MESH_PARTICLES dataStructure;
};


// Property identifier
constexpr int concentration = 1;
constexpr int accumulator = 2;




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







#endif //OPENFPM_PDATA_TP_ALGORITHM_HPP
