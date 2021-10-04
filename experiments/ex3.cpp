//
// Created by landfried on 17.09.21.
//

#include <iostream>
#include <Vector/vector_dist.hpp>


template <typename ParticleType>
class ParticleData {

    float r_cut = 0.3;
    Box<1,float> box;
    size_t bc[1]={PERIODIC};
    Ghost<1,float> ghost;

public:
    vector_dist<1, float, ParticleType> vd;

    ParticleData() : box({0.0},{5.0}), ghost(r_cut), vd(5,box,bc,ghost) {}

    template<unsigned int id> inline auto getProp(vect_dist_key_dx p) -> decltype(vd.template getProp<id>(p)) {
        return vd.template getProp<id>(p);
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

    size_t getID() {
        return key.getKey();
    }

};

template <typename ParticleType>
class ParticleMethod {

    int iteration = 0;

private:
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
    }

    void executeInteraction() {
        auto it2 = particle_data.vd.getDomainIterator();
        while (it2.isNext())
        {
            auto p = it2.get();
            Particle<ParticleType> particle(particle_data, p);

            auto it3 = particle_data.vd.getDomainIterator();
            while (it3.isNext())
            {
                auto n = it3.get();
                Particle<ParticleType> neighbor_particle(particle_data, n);

//                std::cout << n.getKey() << ", " << p.getKey() << std::endl;

                if (p.getKey() != n.getKey()) {
                    // call (overriden) interact method
                    interact(particle, neighbor_particle);
                }

                ++it3;
            }
            ++it2;
        }
    }

protected:
    ParticleData<ParticleType> particle_data;

    virtual void evolve(Particle<ParticleType> particle) {}
    virtual void interact(Particle<ParticleType> particle, Particle<ParticleType> neighbor) {}

public:
    void run() {

        std::cout << "Iteration " << iteration << std::endl;

        executeInteraction();
        executeEvolution();

//        particle_data.vd.deleteGhost();
        particle_data.vd.write_frame("particles",iteration);
//        vd.ghost_get<>();

        iteration++;
    }
};


typedef aggregate<float, float[2]> particle_type;

class TestPM : public ParticleMethod<particle_type> {

    static constexpr int position = 0;
    static constexpr int concentration = 1;

    void evolve(Particle<particle_type> particle) override {
        particle.property<position>() = static_cast<float>(particle.getID());
        particle.property<concentration>()[0] += 2;
        particle.property<concentration>()[1] += 3;
    }

    void interact(Particle<particle_type> particle, Particle<particle_type> neighbor) override {
        std::cout << "interact" << std::endl;
        std::cout << particle.property<position>() << ", " << neighbor.property<position>() << std::endl;

    }

};



int main(int argc, char* argv[]) {

    openfpm_init(&argc,&argv);

    TestPM pm;

    for (int i = 0; i < 5; ++i) {
        pm.run();
    }

    openfpm_finalize();

    return 0;
}
