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
class ParticleMethod {
public:
    typedef ParticleType particleType;
    virtual void evolve(Particle<ParticleType> particle) {}
    virtual void interact(Particle<ParticleType> particle, Particle<ParticleType> neighbor) {}
};

template <typename ParticleMethodType>
class Transition {

protected:

    typedef typename ParticleMethodType::particleType ParticleType;
    ParticleMethodType particleMethod;

//    explicit Transition(ParticleMethod<ParticleType> particleMethod_in) : particleMethod(particleMethod_in), particleData() {}

    int iteration = 0;


    void executeInitialization(ParticleData<ParticleType> &particleData) {
        size_t sz[1] = {10};
        auto it2 = particleData.vd.getGridIterator(sz);
        while (it2.isNext())
        {
            particleData.vd.add();
            auto node = it2.get();
            particleData.vd.getLastPos()[0] = node.get(0) * it2.getSpacing(0);
            particleData.vd.template getLastProp<0>() = node.get(0);
            ++it2;
        }
    }

    void executeEvolution(ParticleData<ParticleType> &particleData) {
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particleData, p);
            // call (overriden) evolve method
            particleMethod.evolve(particle);
//            particleData.vd.getPos(p)[0] = particleData.vd.template getProp<0>(p);
            ++it2;
        }

//        particleData.vd.template ghost_get<>();

    }

    virtual void executeInteraction(ParticleData<ParticleType> &particleData) {
        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particleData, p);

            auto it = particleData.vd.getDomainAndGhostIterator();
            while (it.isNext()) {
                Particle<ParticleType> neighbor(particleData, it.get());
                if (particle != neighbor) {
//                    std::cout << particle.template property<0>() << " neighbor prop 0 " << neighbor.getParticleData().vd.template getProp<0>(neighbor.getID()) << std::endl;
                    particleMethod.interact(particle, neighbor);
                }
                ++it;
            }
            ++it2;
        }
    }

public:

    void initialize(ParticleData<ParticleType> &particleData) {
        executeInitialization(particleData);
//        particleData.vd.map();
//        particleData.vd.template ghost_get<0, 1>();
    }

    void run(ParticleData<ParticleType> &particleData) {
/*
        auto & vcl = create_vcluster();
        if (vcl.getProcessUnitID() == 0) {
            std::cout << "Iteration " << iteration << std::endl;
        }*/

        particleData.vd.map();
        particleData.vd.template ghost_get<0, 1>();

        executeInteraction(particleData);
        executeEvolution(particleData);

        particleData.vd.deleteGhost();
        particleData.vd.write_frame("particles",iteration);

        iteration++;
    }
};

template <typename ParticleMethodType>
class TransitionCellList : public Transition<ParticleMethodType>{
    using typename Transition<ParticleMethodType>::ParticleType;

    CELL_MEMBAL(1, float) cellList;

    void executeInteraction(ParticleData<ParticleType> &particleData) override {
        particleData.vd.template updateCellList(cellList);

        auto it2 = particleData.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particleData, p);

            auto it = cellList.template getNNIterator<NO_CHECK>(cellList.getCell(particleData.vd.getPos(p)));

//            auto it = this->particleData.vd.getDomainAndGhostIterator();
            while (it.isNext()) {
                Particle<ParticleType> neighbor(particleData, it.get());
                if (particle != neighbor) {
//                    std::cout << particle.template property<0>() << " neighbor prop 0 " << neighbor.getParticleData().vd.template getProp<0>(neighbor.getID()) << std::endl;
//                    std::cout << "CellList" << std::endl;
                    this->particleMethod.interact(particle, neighbor);
                }
                ++it;
            }
            ++it2;
        }
    }

public:
    explicit TransitionCellList(ParticleData<ParticleType> &particleData) : Transition<ParticleMethodType>(), cellList(particleData.vd.template getCellList<CELL_MEMBAL(1, float)>(0.5)) {}

};


typedef aggregate<float, float[2]> particle_type;
class TestPM : public ParticleMethod<particle_type> {
public:
    static constexpr int position = 0;
    static constexpr int concentration = 1;

    void evolve(Particle<particle_type> particle) override {

        auto & vcl = create_vcluster();
        std::cout << " evolve cpu " << vcl.getProcessUnitID() << ", particle " << particle.property<position>() << ", position " <<  particle.position()[0]  << std::endl;

//        particle.property<position>() = static_cast<float>(particle.getID());
        particle.property<concentration>()[0] += 2;
        particle.property<concentration>()[1] += 3;
    }

    void interact(Particle<particle_type> particle, Particle<particle_type> neighbor) override {
        auto & vcl = create_vcluster();
        std::cout << " interact cpu " << vcl.getProcessUnitID() << ", particle " << particle.property<position>()
                << ", neighbor " << neighbor.property<position>() << " (" << neighbor.position()[0] << ")" << std::endl;
        particle.property<position>() += 0.1f;

//        std::cout << particle.property<position>() << ", " << neighbor.property<position>() << std::endl;
    }

};



int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    ParticleData<TestPM::particleType> particleData;

    TransitionCellList<TestPM> transition(particleData);
    Transition<TestPM> transition2();

    transition.initialize(particleData);

    for (int i = 0; i < 1; ++i) {
        transition.run(particleData);
    }

    openfpm_finalize();

    return 0;
}
