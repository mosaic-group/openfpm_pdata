//
// Created by landfried on 06.12.21.
//

#ifndef OPENFPM_PDATA_PARTICLEDATA_HPP
#define OPENFPM_PDATA_PARTICLEDATA_HPP

#include <Vector/vector_dist.hpp>
#include "GlobalVar.hpp"
#include "DataContainer.hpp"
template <typename ParticleMethodType, typename SimulationParametersType>
class ParticleData {

    using ParticleSignatureType = typename ParticleMethodType::ParticleSignature;
    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    SimulationParametersType simulationParameters;

    using ContainerFactoryType = DataContainerFactory<ParticleSignatureType>;
    ContainerFactoryType dataContainerFactory;
    using DataContainerType = typename ContainerFactoryType::ContainerType;
    using DataStructureType = typename ContainerFactoryType::ContainerType::DataStructureType;

public:

    DataContainerType dataContainer;

    ParticleData() : dataContainer(dataContainerFactory.template createContainer<SimulationParametersType>(simulationParameters))
    {
        dataContainer.printType();
    }

    DataStructureType& getContainer() {
        return dataContainer.getContainer();
    }

    DataContainerType& getDataContainer() {
        return dataContainer;
    }


};


#endif //OPENFPM_PDATA_PARTICLEDATA_HPP
