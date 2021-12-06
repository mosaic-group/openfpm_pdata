//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLE_HPP
#define OPENFPM_PDATA_PARTICLE_HPP

#include <Vector/vector_dist.hpp>
#include "ParticleData.hpp"

template <typename ParticleType>
class Particle {
    ParticleData<ParticleType>& particle_data;
    vect_dist_key_dx key;

public:
    Particle(ParticleData<ParticleType>& particle_data_in, vect_dist_key_dx key_in) : particle_data(particle_data_in), key(key_in) {}

    template<unsigned int id> inline auto property() -> decltype(particle_data.template getProp<id>(key)) {
        return particle_data.template getProp<id>(key);
    }

    inline auto position() -> decltype(particle_data.getPos(key)) {
        return particle_data.getPos(key);
    }

    size_t getID() {
        return key.getKey();
    }

    vect_dist_key_dx getKey() {
        return key;
    }

    ParticleData<ParticleType>& getParticleData() {
        return particle_data;
    }

    bool operator== (Particle<ParticleType> rhs) {
        return getID() == rhs.getID();
    }

    bool operator!= (Particle<ParticleType> rhs) {
        return getID() != rhs.getID();
    }
};


#endif //OPENFPM_PDATA_PARTICLE_HPP
