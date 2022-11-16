//
// Created by landfried on 14.03.22.
//

#ifndef OPENFPM_PDATA_INSTANCE_HPP
#define OPENFPM_PDATA_INSTANCE_HPP

#include "iostream"
#include "ParticleData.hpp"
#include "../../openfpm_numerics/src/Draw/DrawParticles.hpp"
#include "Particle.hpp"

template <typename ParticleMethodType, typename SimulationParametersType>
class Instance {

protected:

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    typedef Point<dimension, PositionType> PointType;

    ParticleData<ParticleMethodType, SimulationParametersType> &particleData;

    void addParticle() {
        particleData.getOpenFPMContainer().add();
    }

    auto position() {
        return particleData.getOpenFPMContainer().getLastPos();
    }

    template<int dimension>
    auto property() -> decltype(particleData.getOpenFPMContainer().template getLastProp<dimension>()) {
        return particleData.getOpenFPMContainer().template getLastProp<dimension>();
    }

    template <unsigned int dimension, typename PositionType>
    auto boxIterator(Point<dimension, PositionType> pointMin, Point<dimension, PositionType> pointMax, size_t (& sz)[dimension]) {
        Box<dimension, PositionType> drawBox(pointMin, pointMax);
        Box<dimension, PositionType> domain (particleData.simulationParameters.domainMin, particleData.simulationParameters.domainMax);

        auto boxIterator = DrawParticles::DrawBox(particleData.getOpenFPMContainer(),sz,domain,drawBox);
        return boxIterator;
    }

    template <unsigned int dimension, typename PositionType>
    auto skinIterator(Point<dimension, PositionType> pointMin, Point<dimension, PositionType> pointMax, PositionType width,  size_t (& sz)[dimension]) {
        Box<dimension, PositionType> drawOuterBox(pointMin, pointMax);
        Box<dimension, PositionType> drawInnerBox(pointMin + width, pointMax - width);
        Box<dimension, PositionType> domain (particleData.simulationParameters.domainMin, particleData.simulationParameters.domainMax);

        auto boxIterator = DrawParticles::DrawSkin(particleData.getOpenFPMContainer(),sz,domain,drawInnerBox, drawOuterBox);
        return boxIterator;
    }

    template <unsigned int dimension, typename PositionType>
    auto skinOpenIterator(Point<dimension, PositionType> pointMin, Point<dimension, PositionType> pointMax, PositionType width,  size_t (& sz)[dimension]) {
        Point<dimension, PositionType> upperCutoutPoint = pointMax - width;
        upperCutoutPoint[dimension - 1] = pointMax[dimension - 1];

        Box<dimension, PositionType> drawOuterBox(pointMin, pointMax);
        Box<dimension, PositionType> drawInnerBox(pointMin + width, upperCutoutPoint);
        Box<dimension, PositionType> domain (particleData.simulationParameters.domainMin, particleData.simulationParameters.domainMax);

        auto boxIterator = DrawParticles::DrawSkin(particleData.getOpenFPMContainer(),sz,domain,drawInnerBox, drawOuterBox);
        return boxIterator;
    }


    // Utility methods

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()

    float normalDistribution(float mean, float stdev) {
        std::normal_distribution<> dis(mean, stdev);
//        std::uniform_real_distribution<> dis(0, 10);
        return dis(gen);
    }

public:

    Instance(ParticleData<ParticleMethodType, SimulationParametersType> &particleData_in) : particleData(particleData_in) {

    }

    virtual void freePlacement() {}

    virtual void shapePlacement() {}

    //virtual void setGlobalVariable() {}

    virtual void initialization(Particle<ParticleSignatureType> particle) {}


};


#endif //OPENFPM_PDATA_INSTANCE_HPP
