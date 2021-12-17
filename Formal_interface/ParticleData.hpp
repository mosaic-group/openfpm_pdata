//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEDATA_HPP
#define OPENFPM_PDATA_PARTICLEDATA_HPP

#include <Vector/vector_dist.hpp>
#include "ParticleDataReference.hpp"
#include "Domain.hpp"


template <typename ParticleMethodType>
class ParticleData {

    typedef typename ParticleMethodType::particleType ParticleType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    float r_cut = 0.5;
//    Box<ParticleMethodType::spaceDimension,float> box;
//    size_t bc[ParticleMethodType::spaceDimension] = {NON_PERIODIC};
    Ghost<ParticleMethodType::spaceDimension,float> ghost;
    BoundaryCondition<ParticleMethodType::spaceDimension> bc;

public:
    vector_dist<ParticleMethodType::spaceDimension, PositionType , ParticleType> vd;

    ParticleData() : /*box(ParticleMethodType::domainMin,ParticleMethodType::domainMax),*/ ghost(r_cut),
        vd(0, getDomain<dimension, PositionType>(ParticleMethodType::domainMin,ParticleMethodType::domainMax), bc.periodic ,ghost) {

    }

    template<unsigned int id> inline auto getProp(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
    }

    inline auto getPos(vect_dist_key_dx p) -> decltype(vd.getPos(p)) {
        return vd.getPos(p);
    }

    ParticleDataReference<dimension, PositionType, ParticleType> getReference() {
        ParticleDataReference<dimension, PositionType, ParticleType> particleDataReference(&vd);
        return particleDataReference;
    }
};

class ParticleKind {};

template <typename ParticleMethodType, typename = int>
class PD {
public:
    PD() {
        std::cout << "Error particle kind not defined" << std::endl;
    }
};

template <typename ParticleMethodType>
class PD <ParticleMethodType, decltype((void) ParticleMethodType::meshParticles, 0)> {
public:
    PD() {
        std::cout << "PD mesh" << std::endl;
    }
};

template <typename ParticleMethodType>
class PD <ParticleMethodType, decltype((void) ParticleMethodType::freeParticles, 0)> {
public:
    PD() {
        std::cout << "PD free" << std::endl;
    }
};

#endif //OPENFPM_PDATA_PARTICLEDATA_HPP
