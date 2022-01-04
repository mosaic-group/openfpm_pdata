//
// Created by landfried on 21.12.21.
//

#ifndef OPENFPM_PDATA_INITIALCONDITION_HPP
#define OPENFPM_PDATA_INITIALCONDITION_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"

template <typename ParticleMethodType>
class InitialCondition {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

public:
//    constexpr static PositionType domainMin[dimension];
//    constexpr static PositionType domainMax[dimension];


//    virtual void initialization(ParticleRef<ParticleMethodType::dimension, typename ParticleMethodType::PositionType, typename ParticleMethodType::ParticleType> particle) {}

};

#endif //OPENFPM_PDATA_INITIALCONDITION_HPP
