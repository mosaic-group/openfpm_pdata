//
// Created by landfried on 18.01.22.
//

#ifndef OPENFPM_PDATA_DATACONTAINER_HPP
#define OPENFPM_PDATA_DATACONTAINER_HPP

#include <iostream>
#include <Grid/grid_dist_id.hpp>
//#include <Grid/grid_dist_key.hpp>
#include <Vector/vector_dist.hpp>


struct FREE_PARTICLES {};
struct MESH_PARTICLES {};

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


    DataContainer_VectorDist(int numberParticles, Box<dimension, PositionType> domain, const periodicity<dimension> boundaryConditions) :
        ghost(0.5),
        vd(numberParticles, domain, boundaryConditions.bc, ghost) {}

    void printType() override {
        std::cout << "vector_dist" << std::endl;
    }

    template<unsigned int id>
    inline auto property(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
    }

    inline auto position(vect_dist_key_dx p) -> decltype(vd.getPos(p)) {
        return vd.getPos(p);
    }

    vector_dist<dimension, PositionType, PropertyType>& getContainer() {
        return vd;
    }

    void n_prop() {
        constexpr int n = decltype(vd)::value_type::max_prop;
    }

    void deleteGhost() override {
        vd.deleteGhost();
    }

    bool write_frame(std::string out, size_t iteration, int opt = VTK_WRITER) override {
        return vd.write_frame(out, iteration, opt);
    }

};

template <typename ParticleSignatureType>
class DataContainer_GridDist : DataContainer<ParticleSignatureType> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

public:    typedef grid_dist_id<dimension, PositionType, PropertyType> DataStructureType;

private:
    Ghost<2,long int> ghost;

    DataStructureType grid;
//    grid_dist_id<dimension, PositionType , PropertyType> grid;

public:

    DataContainer_GridDist(const size_t (&meshSize)[dimension], Box<dimension, PositionType> domain, const periodicity<dimension> boundaryConditions) :
        ghost(1),
        grid(meshSize,domain,ghost,boundaryConditions) {}

    void printType() override {
        std::cout << "grid_dist" << std::endl;
    }

    template<unsigned int id>
    inline auto property(grid_dist_key_dx<dimension> p) -> decltype(grid.template get<id>(p)) {
        return grid.template get<id>(p);
    }

    inline auto position(grid_dist_key_dx<dimension> p) -> decltype(grid.getPos(p)) {
        return grid.getPos(p);
    }

    grid_dist_id<dimension, PositionType, PropertyType>& getContainer() {
        return grid;
    }

    void deleteGhost() override {}

    bool write_frame(std::string out, size_t iteration, int opt = VTK_WRITER | FORMAT_BINARY) override {
        return grid.write_frame(out,iteration, opt);
    }

};

template<typename ParticleSignatureType, typename DataStructureType = typename ParticleSignatureType::dataStructure>
struct DataContainerFactory {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    DataContainer<ParticleSignatureType> createContainer() {}
};

template<typename ParticleSignatureType>
struct DataContainerFactory<ParticleSignatureType, FREE_PARTICLES> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    typedef DataContainer_VectorDist<ParticleSignatureType> ContainerType;
    typedef vect_dist_key_dx KeyType;


    template<typename SimulationParametersType>
    ContainerType createContainer(SimulationParametersType& simulationParameters) {
        ContainerType newContainer(simulationParameters.numberParticles,
                                   Box<dimension, PositionType>(simulationParameters.domainMin, simulationParameters.domainMax),
                                           simulationParameters.boundaryConditions);
        return newContainer;
    }
};

template<typename ParticleSignatureType>
struct DataContainerFactory<ParticleSignatureType, MESH_PARTICLES> {

    static constexpr int dimension = ParticleSignatureType::dimension;
    using PositionType = typename ParticleSignatureType::position;
    using PropertyType = typename ParticleSignatureType::properties;

    typedef DataContainer_GridDist<ParticleSignatureType> ContainerType;
    typedef grid_dist_key_dx<dimension> KeyType;

    template<typename SimulationParametersType>
    ContainerType createContainer(SimulationParametersType& simulationParameters) {
        ContainerType newContainer(simulationParameters.meshSize, Box<dimension, PositionType>(simulationParameters.domainMin, simulationParameters.domainMax),
                simulationParameters.boundaryConditions);
        return newContainer;
    }
};




#endif //OPENFPM_PDATA_DATACONTAINER_HPP
