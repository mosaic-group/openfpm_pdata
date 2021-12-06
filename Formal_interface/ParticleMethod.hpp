//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEMETHOD_HPP
#define OPENFPM_PDATA_PARTICLEMETHOD_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"


template <typename ParticleType>
class ParticleMethod {
public:
    typedef ParticleType particleType;
    virtual void evolve(Particle<ParticleType> particle) {}
    virtual void interact(Particle<ParticleType> particle, Particle<ParticleType> neighbor) {}
};





#endif //OPENFPM_PDATA_PARTICLEMETHOD_HPP
