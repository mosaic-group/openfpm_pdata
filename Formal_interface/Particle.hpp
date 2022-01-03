//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLE_HPP
#define OPENFPM_PDATA_PARTICLE_HPP

#include <Vector/vector_dist.hpp>
#include "ParticleData.hpp"

template <int dimension, typename PositionType, typename ParticleType>
class Particle {
    vector_dist<dimension, PositionType, ParticleType>& vd_ref;
    vect_dist_key_dx key;

public:
    Particle(vector_dist<dimension, PositionType, ParticleType>& vd_in, vect_dist_key_dx key_in) : vd_ref(vd_in), key(key_in) {}

    template<unsigned int id> inline auto property() -> decltype(vd_ref.template getProp<id>(key)) {
        return vd_ref.template getProp<id>(key);
    }

    inline auto position() -> decltype(vd_ref.getPos(key)) {
        return vd_ref.getPos(key);
    }

    size_t getID() {
        return key.getKey();
    }

    vect_dist_key_dx getKey() {
        return key;
    }

    bool operator== (Particle<dimension, PositionType, ParticleType> rhs) {
        return getID() == rhs.getID();
    }

    bool operator!= (Particle<dimension, PositionType, ParticleType> rhs) {
        return getID() != rhs.getID();
    }
};


#endif //OPENFPM_PDATA_PARTICLE_HPP
