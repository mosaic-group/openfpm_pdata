//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEDATA_HPP
#define OPENFPM_PDATA_PARTICLEDATA_HPP

#include <Vector/vector_dist.hpp>
#include "GlobalVar.hpp"

template <typename ParticleMethodType, typename SimulationParametersType>
class ParticleData {

    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    typedef typename ParticleMethodType::globalVarType GlobalVarType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    float r_cut = 0.5;
    Ghost<dimension, PositionType> ghost;

    SimulationParametersType simulationParameters;


public:
    vector_dist<dimension, PositionType, PropertyType> vd;
    GlobalVar<GlobalVarType> globalVar;

    ParticleData() : ghost(r_cut),
        vd(simulationParameters.numberParticles,
           Box<dimension, PositionType>(simulationParameters.domainMin, simulationParameters.domainMax),
           simulationParameters.boundaryConditions ,ghost) {}

    template<unsigned int id> inline auto getProp(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
    }

    inline auto getPos(vect_dist_key_dx p) -> decltype(vd.getPos(p)) {
        return vd.getPos(p);
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
