//
// Created by landfried on 17.09.21.
//

#include <iostream>
#include <Vector/vector_dist.hpp>

constexpr int NEIGHBOR_ALL = 1;
constexpr int NEIGHBOR_CELLLIST = 2;

template <typename ParticleType>
class ParticleData {

    float r_cut = 0.5;
    Box<1,float> box;
    size_t bc[1]={PERIODIC};
    Ghost<1,float> ghost;

public:
    vector_dist<1, float, ParticleType> vd;

    ParticleData() : box({0.0},{5.0}), ghost(r_cut), vd(0,box,bc,ghost) {}

    template<unsigned int id> inline auto getProp(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
    }

    inline auto getPos(vect_dist_key_dx p) -> decltype(vd.getPos(p)) {
        return vd.getPos(p);
    }
};

template <typename ParticleType>
class Particle {
    ParticleData<ParticleType>& particle_data;
    vect_dist_key_dx key;

public:
    Particle(ParticleData<ParticleType>& particle_data_in, vect_dist_key_dx key_in) : particle_data(particle_data_in), key(key_in) {}

    template<unsigned int id> inline auto property() -> decltype(particle_data.template getProp<id>(key)) {
        return particle_data.template getProp<id>(key);
    }

    inline auto position() -> decltype(particle_data.getPos(key)) {
        return particle_data.getPos(key);
    }

    size_t getID() {
        return key.getKey();
    }

    vect_dist_key_dx getKey() {
        return key;
    }

    ParticleData<ParticleType>& getParticleData() {
        return particle_data;
    }

    bool operator== (Particle<ParticleType> rhs) {
        return getID() == rhs.getID();
    }

    bool operator!= (Particle<ParticleType> rhs) {
        return getID() != rhs.getID();
    }
};

template <typename ParticleType>
class BaseParticleIterator {
protected:
    ParticleData<ParticleType>& particle_data;

public:
    explicit BaseParticleIterator(ParticleData<ParticleType>& particle_data_in) : particle_data(particle_data_in) {}

    virtual bool hasNext() {
        std::cout << "base iterator" << std::endl;
        return false; }
    virtual Particle<ParticleType> getParticle();
    virtual void next() {}
};

template <typename ParticleType>
class ParticleIteratorVectorDist : public BaseParticleIterator<ParticleType> {
private:
    vector_dist_iterator iterator;

public:
    ParticleIteratorVectorDist(ParticleData<ParticleType>& particle_data_in, const vector_dist_iterator& iterator_in) :
        BaseParticleIterator<ParticleType>(particle_data_in), iterator(iterator_in) {}

    bool hasNext() override {
        return iterator.isNext();
    }

    Particle<ParticleType> getParticle() override {
        Particle<ParticleType> particle(this->particle_data, iterator.get());
        return particle;
    }

    void next() override {
        ++iterator;
    }
};

template <typename ParticleType>
class Neighborhood {

protected:
    ParticleData<ParticleType>& particle_data;

public:

    explicit Neighborhood(ParticleData<ParticleType>& particle_data_in) : particle_data(particle_data_in) {}
    virtual void iterationSetup() {}
    virtual void initialize(Particle<ParticleType> particle) {}
    virtual bool isNext() { return false; }
    virtual Particle<ParticleType> getNeighbor();
    virtual void next() {}
    virtual BaseParticleIterator<ParticleType>* getIterator(Particle<ParticleType> particle);

};

template <typename ParticleType>
class NeighborhoodAll : public Neighborhood<ParticleType> {
private:
    vector_dist_iterator iterator;

public:

    explicit NeighborhoodAll(ParticleData<ParticleType>& particle_data_in) : Neighborhood<ParticleType>(particle_data_in), iterator(particle_data_in.vd.getDomainIterator()) {}

    void initialize(Particle<ParticleType> particle) override {
//        iterator = particle.getParticleData().vd.getDomainIterator();
        iterator = particle.getParticleData().vd.getDomainAndGhostIterator();
    }

    bool isNext() override {
        return iterator.isNext();
    }

    Particle<ParticleType> getNeighbor() override {
        Particle<ParticleType> neighbor(this->particle_data, iterator.get());
        return neighbor;
    }

    void next() override {
        ++iterator;
    }

    ParticleIteratorVectorDist<ParticleType>* getIterator(Particle<ParticleType> particle) override {
        return new ParticleIteratorVectorDist<ParticleType>(this->particle_data, this->particle_data.vd.getDomainAndGhostIterator());
    }

};

/*
template <typename ParticleType>
class NeighborhoodCellList : public Neighborhood<ParticleType> {
private:
//    ParticleData<ParticleType> particle_data;
    CELL_MEMBAL(1, float) cellList;
    CellNNIterator<1, CellList<1, float, Mem_bal<unsigned long>, shift<1, float>, openfpm::vector<Point<1, float>, HeapMemory, memory_traits_lin, openfpm::grow_policy_double, 2> >, 3, 1> cliterator;
//    CellNNIterator<1,CellList<1,float,Mem_bal<>,no_transform<1,float>,openfpm::vector<Point<1,float>>>,(int)openfpm::math::pow(3,1), NO_CHECK> cliterator;


public:

    NeighborhoodCellList(ParticleData<ParticleType>& particle_data_in) : Neighborhood<ParticleType>(particle_data_in),
            cellList(particle_data_in.vd.template getCellList<CELL_MEMBAL(1, float)>(0.5)),
            cliterator(cellList.template getNNIterator(cellList.getCell({0}))) {}

    void iterationSetup() override {
        this->particle_data.vd.template updateCellList(cellList);
    }

    void initialize(Particle<ParticleType> particle) override {
        // iterator zurückgeben
        // gemeinsames iterator interface für CellNNIterator und vector_dist_iterator?
        auto it = cellList.template getNNIterator<NO_CHECK>(cellList.getCell({0}));

        std::cout << "type " << boost::core::demangle(typeid(it).name()) << std::endl;
//        cliterator = it;
//        cliterator = cellList.template getNNIterator<NO_CHECK>(cellList.getCell({0}));
//        cliterator = cellList.template getNNIterator(cellList.getCell(this->particle_data.vd.getPos(particle.getKey())));
    }

    bool isNext() override {
        return cliterator.isNext();
    }

    Particle<ParticleType> getNeighbor() override {
        Particle<ParticleType> neighbor(this->particle_data, cliterator.get());
        return neighbor;
    }

    void next() override {
        ++cliterator;
    }

};
*/

template <typename ParticleType, typename NeighborFunction = NeighborhoodAll<ParticleType>, int Neighborhood = NEIGHBOR_ALL>
class ParticleMethod {


private:

    void executeInitialization() {

        size_t sz[1] = {10};
        auto it2 = particle_data.vd.getGridIterator(sz);
        while (it2.isNext())
        {
            particle_data.vd.add();
            auto node = it2.get();
            particle_data.vd.getLastPos()[0] = node.get(0) * it2.getSpacing(0);
            particle_data.vd.template getLastProp<0>() = node.get(0);
//            Particle<ParticleType> particle(particle_data, p);
            // call (overriden) initialization method
//            initialization(particle);
            ++it2;
        }

        particle_data.vd.map();
        particle_data.vd.template ghost_get<>();

    }

    void executeEvolution() {
        auto it2 = particle_data.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particle_data, p);
            // call (overriden) evolve method
            evolve(particle);
//            particle_data.vd.getPos(p)[0] = particle_data.vd.template getProp<0>(p);
            ++it2;
        }

//        particle_data.vd.template ghost_get<>();

    }

    void executeInteraction() {
        auto it2 = particle_data.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particle_data, p);

            auto it = particle_data.vd.getDomainAndGhostIterator();
            while (it.isNext()) {
                Particle<ParticleType> neighbor(particle_data, it.get());
                if (particle != neighbor) {
//                    std::cout << particle.template property<0>() << " neighbor prop 0 " << neighbor.getParticleData().vd.template getProp<0>(neighbor.getID()) << std::endl;
                    interact(particle, neighbor);
                }
                ++it;
            }
/*
            neighborhood.initialize(particle);
            while (neighborhood.isNext()) {
                Particle<ParticleType> neighbor = neighborhood.getNeighbor();
                if (particle != neighbor) {
//                    std::cout << particle.template property<0>() << " neighbor prop 0 " << neighbor.getParticleData().vd.template getProp<0>(neighbor.getID()) << std::endl;
                    interact(particle, neighbor);
                }
                neighborhood.next();
            }
*/

            ++it2;
        }
    }

protected:

    ParticleMethod() : neighborhood(particle_data) {}

    ParticleData<ParticleType> particle_data;
    NeighborFunction neighborhood;

    int iteration = 0;
    virtual void initialization(Particle<ParticleType> particle) {}
    virtual void evolve(Particle<ParticleType> particle) {}

    virtual void interact(Particle<ParticleType> particle, Particle<ParticleType> neighbor) {}

public:

    void initialize() {
        executeInitialization();
        particle_data.vd.template ghost_get<0, 1>();
    }

    void run() {
/*
        auto & vcl = create_vcluster();
        if (vcl.getProcessUnitID() == 0) {
            std::cout << "Iteration " << iteration << std::endl;
        }*/

        executeInteraction();
        executeEvolution();

        particle_data.vd.deleteGhost();
        particle_data.vd.write_frame("particles",iteration);
        particle_data.vd.template ghost_get<0, 1>();

        iteration++;
    }
};


typedef aggregate<float, float[2]> particle_type;
class TestPM : public ParticleMethod<particle_type> {

    static constexpr int position = 0;
    static constexpr int concentration = 1;

    void evolve(Particle<particle_type> particle) override {

        auto & vcl = create_vcluster();
        std::cout << this->iteration << " evolve cpu " << vcl.getProcessUnitID() << ", particle " << particle.property<position>() << ", position " <<  particle.position()[0]  << std::endl;

//        particle.property<position>() = static_cast<float>(particle.getID());
        particle.property<concentration>()[0] += 2;
        particle.property<concentration>()[1] += 3;
    }

    void interact(Particle<particle_type> particle, Particle<particle_type> neighbor) override {
        auto & vcl = create_vcluster();
        std::cout << this->iteration << " interact cpu " << vcl.getProcessUnitID() << ", particle " << particle.property<position>()
                << ", neighbor " << neighbor.property<position>() << " (" << neighbor.position()[0] << ")" << std::endl;

//        std::cout << particle.property<position>() << ", " << neighbor.property<position>() << std::endl;
    }

};



int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    TestPM pm;
    pm.initialize();

    for (int i = 0; i < 1; ++i) {
        pm.run();
    }

    openfpm_finalize();

    return 0;
}
