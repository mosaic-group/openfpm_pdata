//
// Created by landfried on 18.01.22.
//

#ifndef OPENFPM_PDATA_DATACONTAINER_HPP
#define OPENFPM_PDATA_DATACONTAINER_HPP

#include <iostream>

struct FREE_PARTICLES {};
struct MESH_PARTICLES {};

template <int dimension, typename PositionType, typename PropertyType>
class DataContainer {

public:
    virtual void printType() = 0;

};

template <int dimension, typename PositionType, typename PropertyType>
class DataContainer_VectorDist : DataContainer<dimension, PositionType, PropertyType> {

    vector_dist<dimension, PositionType, PropertyType> vd;

public:

    DataContainer_VectorDist(int numberParticles, Box<dimension, PositionType> domain, const size_t (&boundaryConditions)[dimension], Ghost<dimension, PositionType> ghost) :
        vd(numberParticles, domain, boundaryConditions, ghost) {}

    void printType() override {
        std::cout << "vector_dist" << std::endl;
    }

    template<unsigned int id> inline auto property(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
    }

    inline auto position(vect_dist_key_dx p) -> decltype(vd.getPos(p)) {
        return vd.getPos(p);
    }

};

template <int dimension, typename PositionType, typename PropertyType>
class DataContainer_GridDist : DataContainer<dimension, PositionType, PropertyType> {
public:
    void printType() override {
        std::cout << "grid_dist" << std::endl;
    }

};





template<typename DataStructureType, typename ParticleMethodType, typename SimulationParametersType>
struct DataContainerFactory {
    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    DataContainer<dimension, PositionType, PropertyType> getContainer() {}
};

template<typename ParticleMethodType, typename SimulationParametersType>
struct DataContainerFactory<FREE_PARTICLES, ParticleMethodType, SimulationParametersType> {
    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    typedef DataContainer_VectorDist<dimension, PositionType, PropertyType> ContainerType;

    DataContainer_VectorDist<dimension, PositionType, PropertyType> getContainer(SimulationParametersType& simulationParameters) {
        Ghost<dimension, PositionType> ghost(0.5);
        DataContainer_VectorDist<dimension, PositionType, PropertyType> newContainer(simulationParameters.numberParticles,
                                                                                     Box<dimension, PositionType>(simulationParameters.domainMin, simulationParameters.domainMax),
                                                                                     simulationParameters.boundaryConditions, ghost);
        return newContainer;
    }
};

template<typename ParticleMethodType, typename SimulationParametersType>
struct DataContainerFactory<MESH_PARTICLES, ParticleMethodType, SimulationParametersType> {
    typedef typename ParticleMethodType::propertyType PropertyType;
    typedef typename ParticleMethodType::positionType PositionType;
    static constexpr int dimension = ParticleMethodType::spaceDimension;

    typedef DataContainer_GridDist<dimension, PositionType, PropertyType> ContainerType;
    DataContainer_GridDist<dimension, PositionType, PropertyType> getContainer() {
        DataContainer_GridDist<dimension, PositionType, PropertyType> newContainer;
        return newContainer;
    }
};




#endif //OPENFPM_PDATA_DATACONTAINER_HPP
