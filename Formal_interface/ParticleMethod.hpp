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


    ParticleMethod() : gen(rd()) {}

    virtual void evolve(/*GlobalVar<GlobalVarType> globalVar,*/ Particle<dimension, PositionType, PropertyType> particle) {}
    virtual void interact(Particle<dimension, PositionType, PropertyType> particle, Particle<dimension, PositionType, PropertyType> neighbor) {}
    virtual void initialization(Particle<dimension, PositionType, PropertyType> particle) {}

    virtual void evolveGlobalVar() {}
    virtual bool stop() {
        return true;
    }




    // Random number generator

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()

    float normalDistribution(float mean, float stdev) {
        std::normal_distribution<> dis(mean, stdev);
        return dis(gen);
    }

    /*template<int d, typename T>
    Point<d, T> normalDistribution(float mean, float stdev) {
        std::normal_distribution<> dis(mean, stdev);

        Point<d, T> randomPoint();
        for (int i = 0; i < d; ++i) {
            randomPoint()[i] = dis(gen);
        }
        return randomPoint();
    }*/
};






#endif //OPENFPM_PDATA_PARTICLEMETHOD_HPP
