//
// Created by landfried on 21.09.21.
//



template <int propId>
struct Property {
    static const int id = propId;
    constexpr int operator() () { return id; }
};


class Particle_old {
public:
    static const int dimension = 1;
    Property<0> position;
    Property<1> concentration;
};

