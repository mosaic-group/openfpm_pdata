//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEMETHOD_HPP
#define OPENFPM_PDATA_PARTICLEMETHOD_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"
#include "GlobalVar.hpp"

template <int dimension, typename PositionType, typename PropertyType, typename GlobalVarType>
class ParticleMethod {
public:
    typedef PositionType positionType;
    typedef PropertyType propertyType;
    typedef GlobalVarType globalVarType;
    constexpr static int spaceDimension = dimension;
//    constexpr static float domainMin[spaceDimension] = {0.0};
//    constexpr static float domainMax[spaceDimension] = {1.0};
//    constexpr static size_t boundaryCondition[spaceDimension] = {PERIODIC};


    virtual void evolve(/*GlobalVar<GlobalVarType> globalVar,*/ Particle<dimension, PositionType, PropertyType> particle) {}
    virtual void interact(Particle<dimension, PositionType, PropertyType> particle, Particle<dimension, PositionType, PropertyType> neighbor) {}
    virtual void evolveGlobalVar() {}
    virtual bool stop() {
        return true;
    }
};






#endif //OPENFPM_PDATA_PARTICLEMETHOD_HPP
