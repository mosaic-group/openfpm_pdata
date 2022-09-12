//
// Created by landfried on 02.03.22.
//
// source /home/peter/openfpm_vars
// make run 

#ifndef OPENFPM_PDATA_SPH_ALGORITHM_HPP
#define OPENFPM_PDATA_SPH_ALGORITHM_HPP


#include "Vector/vector_dist.hpp"
#include "Formal_interface/Particle.hpp"
#include "Formal_interface/ParticleData.hpp"
#include "Formal_interface/ParticleMethod.hpp"
#include "Formal_interface/Transition.hpp"
#include "Formal_interface/SimulationParameters.hpp"
#include "Formal_interface/InitialCondition.hpp"
#include "Formal_interface/Neighborhood.hpp"
#include "Formal_interface/Instance.hpp"




constexpr int DIMENSION = 3;
typedef double POSITIONTYPE;

// Property identifier
constexpr int positionOld = 0;
constexpr int velocity = 1;
constexpr int velocityOld = 2;
constexpr int density = 3;
constexpr int densityOld = 4;
constexpr int deltaVelocity = 5;
constexpr int deltaDensity = 6;
constexpr int boundary = 7;


struct SPH_ParticleSignature {
    static constexpr int dimension = DIMENSION;
    typedef POSITIONTYPE position;
    typedef aggregate<double[dimension], double[dimension], double[dimension], double, double, double[dimension],double, bool> properties;
    typedef FREE_PARTICLES dataStructure;
};


struct GlobalVariable {

    static double t;
    static double dt;
    static double endT;

    static double particleSpacing;
    static double particleSpacingWater;
    static double mass;
    static Point<DIMENSION, POSITIONTYPE> gravity;
    static double c0;
    static double density0;
    static int gamma;
    static double nu;
    static double h;//characteristic length

    static double phase;
    static int support;
    static double rc;//cutof radius
    static double epsilon;

    static double domain_min[3];
    static double domain_max[3];

    // calculate number of particles in each dimension
    static size_t sz[3];

} g;


template <typename ParticleSignature>
class SPH_ParticleMethod : public ParticleMethod<ParticleSignature> {
    static constexpr int dimension = ParticleSignature::dimension;
    using PositionType = typename ParticleSignature::position;


    double pressure_density2(double density){
        return g.c0*g.c0*g.density0/g.gamma * (pow(density/g.density0, g.gamma)-1)/density/density;
    }
public:


    void interact(Particle<ParticleSignature> particle, Particle<ParticleSignature> neighbor) override {

        Point<dimension, PositionType> n_pos = neighbor.position();
        Point<dimension, PositionType> p_pos = particle.position();
        Point<dimension, PositionType> r_pq = n_pos - p_pos;
        Point<dimension, PositionType> v_pq =  neighbor.template property_vec<velocity>()
                                               -particle.template property_vec<velocity>();
//        double dist2_pq = r_pq.abs2();
        double dist2_pq = p_pos.distance2(n_pos);
        double dist_pq = sqrt(dist2_pq);

        /*
        if (dist_pq > g.rc)
            return;
*/

        double f_pq=  pow(1.0-dist_pq/2.0/g.h , 3);
        double vr = 0;
        vr = scalarProduct(r_pq, v_pq);

        double p_pressure_density2=pressure_density2(particle.template property_vec<density>());
        double q_pressure_density2=pressure_density2(neighbor.template property_vec<density>());
        double interim01= p_pressure_density2+q_pressure_density2 ;
        double interim02 = 10*g.nu/dist2_pq * vr;

        particle.template property_vec<deltaVelocity>() += (interim01-interim02/particle.template property_vec<density>())* r_pq * f_pq;
        neighbor.template property_vec<deltaVelocity>() -= (interim01-interim02/neighbor.template property_vec<density>())* r_pq * f_pq;



double densityChange= vr * f_pq;

particle.template property_vec<deltaDensity>() += densityChange;
neighbor.template property_vec<deltaDensity>() += densityChange;


}

void evolve(Particle<ParticleSignature> particle) override {
double prefact = g.mass* -5.0*21.0/16.0/M_PI/pow(g.h,5);

Point<dimension, PositionType> acceleration = g.gravity + particle.template property_vec<deltaVelocity>()*prefact;

double densityacceleration = particle.template property_vec<deltaDensity>()*prefact;

        if (g.phase==0){
            if (particle.template property_vec<boundary>()==false){
                // fluid

                // move particles half step
                particle.template property_vec<positionOld>()=particle.position_vec();
                particle.position_vec() += g.dt/2.0f*particle.template property_vec<velocity>();

                // change velocity half step
                particle.template property_vec<velocityOld>()=particle.template property_vec<velocity>();
                particle.template property_vec<velocity>()+= g.dt/2.0*acceleration;
            }

            // fluid + boundary

            // change density
            particle.template property_vec<densityOld>()=particle.template property_vec<density>();
            particle.template property_vec<density>()+=g.dt/2.0*densityacceleration;
        }
        else{
            if (particle.template property_vec<boundary>()==false){
                // fluid

                // move particle from original position
                Point<dimension, PositionType> step_acc_half = g.dt  / 2.0 * acceleration;
                Point<dimension, PositionType> step_vel = g.dt * (particle.template property_vec<velocityOld>() + step_acc_half);
                Point<dimension, PositionType> new_pos = particle.template property_vec<positionOld>() + step_vel;
                particle.position_vec() = new_pos;

                // change velocity
                Point<dimension, PositionType> step_acc_full = g.dt * acceleration;
                particle.template property_vec<velocity>() = particle.template property_vec<velocityOld>() + step_acc_full;


//                particle.position_vec() = particle.template property_vec<positionOld>() + g.dt * (particle.template property_vec<velocityOld>()  +  (g.dt  / 2.0) * acceleration);
//                particle.template property_vec<velocity>() = particle.template property_vec<velocityOld>() +g.dt*acceleration;


            }

            // fluid + boundary

            // change density
            particle.template property_vec<density>()= particle.template property_vec<densityOld>() + g.dt*densityacceleration;
        }

        //set to 0 to have a fresh accumulators for the next time step
        particle.template property_vec<deltaVelocity>()=Point<dimension, PositionType> (0.0);
        particle.template property_vec<deltaDensity>()=0.0;
    }

    void evolveGlobalVariable() override {

        if (g.phase==0){
            g.phase=1;
        }
        else{
            g.t+=g.dt;//advance the current time by the time step
            g.phase=0;
        }

        std::cout << "\r" << "t = " << g.t << " (" << int(g.t / g.endT * 100) << "%) " << std::flush;

    }

    bool stop() override {

        // Check simulation time
        if (g.t > g.endT)
            return true;

        return false;
    }
};




#endif //OPENFPM_PDATA_SPH_ALGORITHM_HPP
