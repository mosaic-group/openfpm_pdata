//
// Created by landfried on 21.01.22.
//

#ifndef OPENFPM_PDATA_PARTICLESIGNATURE_HPP
#define OPENFPM_PDATA_PARTICLESIGNATURE_HPP

#include "Vector/vector_dist.hpp"
#include "DataContainer.hpp"


/**
 * Defines the structure of the particle model used in the simulation.
 */
struct ParticleSignature {

    // Number of dimensions of the space
    static constexpr int dimension = 2;

    // Data-type of the space / particle position
    typedef float position;

    // Properties of the particle in their respective order
    typedef aggregate<float[dimension], float[dimension]> properties;

    // Underlying data structure of the particles.
    // FREE_PARTICLES for Lagrangian particles that can change their position
    // MESH_PARTICLES for Eularian particles with fix positions
    typedef FREE_PARTICLES dataStructure;
};

#endif //OPENFPM_PDATA_PARTICLESIGNATURE_HPP
