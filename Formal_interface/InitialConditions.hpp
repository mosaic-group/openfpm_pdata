//
// Created by landfried on 21.12.21.
//

#ifndef OPENFPM_PDATA_INITIALCONDITIONS_HPP
#define OPENFPM_PDATA_INITIALCONDITIONS_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"

template <typename ParticleMethodType>
class InitialConditions {
//    domain size
//    particle type (mesh/free)

    virtual void initialization(ParticleRef<ParticleMethodType::dimension, typename ParticleMethodType::PositionType, typename ParticleMethodType::ParticleType> particle) {}

};

#endif //OPENFPM_PDATA_INITIALCONDITIONS_HPP
