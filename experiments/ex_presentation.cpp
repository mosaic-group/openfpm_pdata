//
// Created by landfried on 17.09.21.
//

#include <iostream>
#include <Vector/vector_dist.hpp>

typedef aggregate<float, float> particle_type;

class DEM : public ParticleMethod<particle_type> {
public:
    static constexpr int position = 0;
    static constexpr int velocity = 1;

    void evolve(Particle<particle_type> particle) override {
        particle.property<position>() += particle.property<velocity>();
    }

    void interact(Particle<particle_type> particle, Particle<particle_type> neighbor) override {
        float particle_velocity = particle.property<velocity>();
        particle.property<velocity>() = neighbor.property<velocity>();
        neighbor.property<velocity>() = particle_velocity;
    }

};

int main(int argc, char* argv[]) {

    ParticleData<DEM::particleType> particleData;
    TransitionCellList<DEM> transition(particleData);

    for (int i = 0; i < 10; ++i) {
        transition.run(particleData);
    }
}
