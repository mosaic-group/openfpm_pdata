//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEMETHOD_HPP
#define OPENFPM_PDATA_PARTICLEMETHOD_HPP

#include <Vector/vector_dist.hpp>
#include "Particle.hpp"
#include "GlobalVar.hpp"

template <typename ParticleSignatureType>
class ParticleMethod {
public:
//    typedef PositionType positionType;
//    typedef PropertyType propertyType;
//    typedef GlobalVarType globalVarType;
//    constexpr static int spaceDimension = dimension;

    typedef ParticleSignatureType ParticleSignature;

    ParticleMethod() : gen(rd()) {}

    virtual void evolve(Particle<ParticleSignatureType> particle) {}
    virtual void interact(Particle<ParticleSignatureType> particle, Particle<ParticleSignatureType> neighbor) {}
//    virtual void initialization(Particle<ParticleSignature> particle) {}

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
