//
// Created by landfried on 21.01.22.
//

#ifndef OPENFPM_PDATA_PARTICLESIGNATURE_HPP
#define OPENFPM_PDATA_PARTICLESIGNATURE_HPP

#include "Vector/vector_dist.hpp"
#include "DataContainer.hpp"

struct ParticleSignature {
    static constexpr int dimension = 2;
    typedef float position;
    typedef aggregate<float[dimension], float[dimension]> properties;
    typedef FREE_PARTICLES dataStructure;
};

#endif //OPENFPM_PDATA_PARTICLESIGNATURE_HPP
