/*
 * vector_dist_cl_performance_tests.hpp
 *
 *
 *  Created on: Mar 22, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_

#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Vector/vector_dist_unit_test.hpp"
#include "Plot/GoogleChart.hpp"
#include <functional>

BOOST_AUTO_TEST_SUITE( vector_dist_cl_perf_test )

///////////////////// INPUT DATA //////////////////////

// Cut-off radiuses. Can be put different number of values
openfpm::vector<float> r_cutoff {0.007,0.02};
// Orders of a curve. Can be put different number of values
openfpm::vector<size_t> orders = {1,2,3,5,7};
// Number of steps of moving the particles
size_t n_moving = 8;
// Moving distance (step size)
double dist = 0.1;
// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The minimal amount of particles
size_t k_min = 15000;
// Ghost part of distributed vector
double ghost_part = 0.01;

///////////////////////////////////////////////////////

// Numbers of particles vector
openfpm::vector<size_t> n_particles;
// Vectors to store the data for 2D
openfpm::vector<openfpm::vector<double>> time_rand;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_hilb;
openfpm::vector<openfpm::vector<double>> time_total_rand;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_total_hilb;
openfpm::vector<openfpm::vector<openfpm::vector<openfpm::vector<double>>>> time_hilb_moved;
// Vectors to store the data for 3D
openfpm::vector<openfpm::vector<double>> time_rand_2;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_hilb_2;
openfpm::vector<openfpm::vector<double>> time_total_rand_2;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_total_hilb_2;
openfpm::vector<openfpm::vector<openfpm::vector<openfpm::vector<double>>>> time_hilb_moved_2;


BOOST_AUTO_TEST_CASE( vector_dist_cl_random_test )
{
	//Benchmark test for 2D and 3D
	vd_cl_random_benchmark<2>(k_start,k_min,ghost_part,r_cutoff,n_particles,time_rand,time_total_rand);
	vd_cl_random_benchmark<3>(k_start,k_min,ghost_part,r_cutoff,n_particles,time_rand_2,time_total_rand_2);
}

BOOST_AUTO_TEST_CASE( vector_dist_cl_hilbert_test )
{
	//Benchmark test for 2D and 3D
	vd_cl_hilbert_benchmark<2>(k_start,k_min,ghost_part,n_moving,dist,r_cutoff,n_particles,orders,time_hilb,time_total_hilb,time_hilb_moved);
	vd_cl_hilbert_benchmark<3>(k_start,k_min,ghost_part,n_moving,dist,r_cutoff,n_particles,orders,time_hilb_2,time_total_hilb_2,time_hilb_moved_2);
}

BOOST_AUTO_TEST_CASE(vector_dist_cl_performance_write_report)
{
	//Write report for 2D and 3D
	vd_cl_performance_write_report<2>(n_moving,r_cutoff,n_particles,orders,time_hilb,time_rand,time_total_hilb,time_total_rand,time_hilb_moved);
	vd_cl_performance_write_report<3>(n_moving,r_cutoff,n_particles,orders,time_hilb_2,time_rand_2,time_total_hilb_2,time_total_rand_2,time_hilb_moved_2);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_ */
