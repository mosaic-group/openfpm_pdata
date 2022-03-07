//
// Created by landfried on 18.01.22.
//

#ifndef OPENFPM_PDATA_DATACONTAINER_HPP
#define OPENFPM_PDATA_DATACONTAINER_HPP

#include <iostream>
#include "Grid/grid_dist_id.hpp"
//#include <Grid/grid_dist_key.hpp>
#include "Vector/vector_dist.hpp"
#include "OperationProxy.hpp"

struct FREE_PARTICLES {};
struct MESH_PARTICLES {};

/**
 * Primary template
 * This class contains the actual OpenFPM data structure.
 * There are 2 child classes for either free particles (using openfpm::vector_dist) or mesh particles (using openfpm::grid_dist_id).
 * The class contains virtual functions that are implemented by the children with respect to the data structure.
 * @tparam ParticleSignatureType
 */
template <typename ParticleSignatureType>
class DataContainer {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:

    virtual void printType() = 0;

    // OpenFPM functions
    virtual void deleteGhost() = 0;
    virtual bool write_frame(std::string out, size_t iteration, int opt) = 0;

};

/**
 * Implementation of DataContainer for free particles.
 * Data is stored in an openfpm::vector_dist.
 * @tparam ParticleSignatureType
 */
template <typename ParticleSignatureType>
class DataContainer_VectorDist : DataContainer<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:    typedef vector_dist<dimension, PositionType, PropertyType> DataStructureType;

private:
    Ghost<dimension, PositionType> ghost;

    DataStructureType vd;

public:

    DataContainer_VectorDist(int numberParticles, Box<dimension, PositionType> domain, const size_t (&boundaryConditions)[dimension], int ghostSize) :
        ghost(ghostSize),
        vd(numberParticles, domain, boundaryConditions, ghost) {}

    void printType() override {
        std::cout << "vector_dist" << std::endl;
    }

    /**
     * Returns a reference to a particle property in the vector_dist data structure.
     * @tparam id Component of the property
     * @param p Referenced particle
     * @return Reference to property field
     */
    template<unsigned int id>
    inline auto property(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
    }

    template<unsigned int id>
    inline auto property_vec(vect_dist_key_dx p) -> OperationProxy<typename std::remove_reference<decltype(property<id>(p))>::type>  {
        OperationProxy<typename std::remove_reference<decltype(property<id>(p))>::type> operationProxy(property<id>(p));
        return operationProxy;
    }

    /**
     * Returns a reference to the position property of a particle in the vector_dist data structure.
     * @param p Referenced particle
     * @return Reference to the position property
     */
    inline auto position(vect_dist_key_dx p) -> decltype(vd.getPos(p)) {
        return vd.getPos(p);
    }

    virtual inline auto position_vec(vect_dist_key_dx p) -> OperationProxy<typename std::remove_reference<decltype(position(p))>::type> {
        OperationProxy<typename std::remove_reference<decltype(position(p))>::type> operationProxy(position(p));
        return operationProxy;
    }


    /**
     * Returns reference to the underlying vector_dist data structure
     * @return vector_dist reference
     */
    vector_dist<dimension, PositionType, PropertyType>& getContainer() {
        return vd;
    }

    /**
     * Returns the number of properties
     */
    void n_prop() {
        constexpr int n = decltype(vd)::value_type::max_prop;
    }

    /**
     * Implementation of deleting the ghost
     */
    void deleteGhost() override {
        vd.deleteGhost();
    }

    /**
     * Writes the current state to a file
     * @param out File name
     * @param iteration Time step (is added to file name)
     * @param opt Output format
     * @return Success of writing to the file
     */
    bool write_frame(std::string out, size_t iteration, int opt = VTK_WRITER) override {
        return vd.write_frame(out, iteration, opt);
    }

};

/**
 * Implementation of DataContainer for mesh particles.
 * Data is stored in an openfpm::grid_dist_id.
 * @tparam ParticleSignatureType
 */
template <typename ParticleSignatureType>
class DataContainer_GridDist : DataContainer<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:    typedef grid_dist_id<dimension, PositionType, PropertyType> DataStructureType;

private:
    Ghost<2,long int> ghost;

    DataStructureType grid;

public:

    DataContainer_GridDist(const size_t (&meshSize)[dimension], Box<dimension, PositionType> domain, const periodicity<dimension> boundaryConditions, int ghostSize) :
        ghost(ghostSize),
        grid(meshSize,domain,ghost,boundaryConditions) {}

    void printType() override {
        std::cout << "grid_dist" << std::endl;
    }

    /**
     * Returns a reference to a particle property in the grid_dist_id data structure.
     * @tparam id Component of the property
     * @param p Referenced particle
     * @return Reference to property field
     */
    template<unsigned int id>
    inline auto property(grid_dist_key_dx<dimension> p) -> decltype(grid.template get<id>(p)) {
        return grid.template get<id>(p);
    }


    template<unsigned int id>
    inline auto property_vec(grid_dist_key_dx<dimension> p) -> OperationProxy<typename std::remove_reference<decltype(property<id>(p))>::type>  {
        OperationProxy<typename std::remove_reference<decltype(property<id>(p))>::type> operationProxy(property<id>(p));
        return operationProxy;
    }


    /**
     * Returns a reference to the position property of a particle in the vector_dist data structure.
     * @param p Referenced particle
     * @return Reference to the position property
     */
    inline auto position(grid_dist_key_dx<dimension> p) -> decltype(grid.getPos(p)) {
        return grid.getPos(p);
    }

    inline auto position_vec(grid_dist_key_dx<dimension> p) {
        return position(p);
    }

    /**
     * Returns reference to the underlying grid_dist_id data structure
     * @return grid_dist_id reference
     */
    grid_dist_id<dimension, PositionType, PropertyType>& getContainer() {
        return grid;
    }

    /**
     * Implementation of deleting the ghost
     */
    void deleteGhost() override {}

    /**
     * Writes the current state to a file
     * @param out File name
     * @param iteration Time step (is added to file name)
     * @param opt Output format
     * @return Success of writing to the file
     */
    bool write_frame(std::string out, size_t iteration, int opt = VTK_WRITER | FORMAT_ASCII) override {
        return grid.write_frame(out,iteration, opt);
    }

};

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
                                           simulationParameters.boundaryConditions, simulationParameters.cutoff_radius);
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




#endif //OPENFPM_PDATA_DATACONTAINER_HPP
