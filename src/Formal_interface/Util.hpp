//
// Created by landfried on 14.11.22.
//

#ifndef OPENFPM_PDATA_UTIL_HPP
#define OPENFPM_PDATA_UTIL_HPP

#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"


template<typename positionType, unsigned int dimension>
positionType abs2(Point<dimension, positionType> point) {
    return point.distance2(point.zero_p());
}

template<typename ParticleSignatureType>
typename ParticleSignatureType::position distance2(Particle<ParticleSignatureType> particle, Particle<ParticleSignatureType> neighbor) {
    Point<ParticleSignatureType::dimension, typename ParticleSignatureType::position> p_pos = particle.position_raw();
    Point<ParticleSignatureType::dimension, typename ParticleSignatureType::position> n_pos = neighbor.position_raw();
    return p_pos.distance2(n_pos);
}

template<typename ParticleSignatureType>
typename ParticleSignatureType::position distance(Particle<ParticleSignatureType> particle, Particle<ParticleSignatureType> neighbor) {
    Point<ParticleSignatureType::dimension, typename ParticleSignatureType::position> p_pos = particle.position_raw();
    Point<ParticleSignatureType::dimension, typename ParticleSignatureType::position> n_pos = neighbor.position_raw();
    return p_pos.distance(n_pos);
}

#endif //OPENFPM_PDATA_UTIL_HPP
