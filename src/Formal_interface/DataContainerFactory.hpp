//
// Created by landfried on 20.09.22.
//

#ifndef OPENFPM_PDATA_DATACONTAINERFACTORY_HPP
#define OPENFPM_PDATA_DATACONTAINERFACTORY_HPP

#include "DataContainer.hpp"


/**
 * Primary template
 * Factory for creating a DataContainer object according to the ParticleSignature.dataStructure type
 * @tparam ParticleSignatureType
 * @tparam DataStructureType Determines which data structure is used to store the particles
 */
template<typename ParticleSignatureType, typename DataStructureType = typename ParticleSignatureType::dataStructure>
struct DataContainerFactory {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    DataContainer<ParticleSignatureType> createContainer() {}
};

/**
 * Factory implementation for creating a DataContainer object for free particles.
 * @tparam ParticleSignatureType
 */
template<typename ParticleSignatureType>
struct DataContainerFactory<ParticleSignatureType, FREE_PARTICLES> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    typedef DataContainer_VectorDist<ParticleSignatureType> ContainerType;
    typedef vect_dist_key_dx KeyType;

    /**
     * Calls constructor with appropriate simulation parameters
     * @tparam SimulationParametersType
     * @param simulationParameters contains parameters used to construct the vector_dist
     * @return DataContainer_VectorDist for free particles
     */
    template<typename SimulationParametersType>
    ContainerType createContainer(SimulationParametersType& simulationParameters) {
        ContainerType newContainer(simulationParameters.numberParticles,
                                   Box<dimension, PositionType>(simulationParameters.domainMin, simulationParameters.domainMax),
                                   simulationParameters.boundaryConditions, simulationParameters.cutoff_radius, simulationParameters.decompositionGranularity);
        return newContainer;
    }
};

/**
 * Factory implementation for creating a DataContainer object for mesh particles
 * @tparam ParticleSignatureType
 */
template<typename ParticleSignatureType>
struct DataContainerFactory<ParticleSignatureType, MESH_PARTICLES> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    typedef DataContainer_GridDist<ParticleSignatureType> ContainerType;
    typedef grid_dist_key_dx<dimension> KeyType;

    /**
     * Calls constructor with appropriate simulation parameters
     * @tparam SimulationParametersType
     * @param simulationParameters contains parameters used to construct the grid_dist_id
     * @return DataContainer_GridDist for mesh particles
     */
    template<typename SimulationParametersType>
    ContainerType createContainer(SimulationParametersType& simulationParameters) {

        // calculate cutoff radius in mesh nodes
        int neighborDistance = std::ceil(simulationParameters.cutoff_radius / simulationParameters.meshSpacing);

        // create grid_dist_id boundary conditions
        periodicity<dimension> boundaryConditions;
        for (int i = 0; i < dimension; ++i)
            boundaryConditions.bc[i] = simulationParameters.boundaryConditions[i];


        ContainerType newContainer(simulationParameters.meshSize, Box<dimension, PositionType>(simulationParameters.domainMin, simulationParameters.domainMax),
                                   boundaryConditions, neighborDistance);
        return newContainer;
    }
};




#endif //OPENFPM_PDATA_DATACONTAINERFACTORY_HPP
