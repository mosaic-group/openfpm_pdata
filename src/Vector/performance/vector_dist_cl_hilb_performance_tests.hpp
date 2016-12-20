/*
 * vector_dist_cl_hilb_performance_tests.hpp
 *
 *  Created on: May 24, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_CL_HILB_PERFORMANCE_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_CL_HILB_PERFORMANCE_TESTS_HPP_

#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "vector_dist_performance_util.hpp"

BOOST_AUTO_TEST_SUITE( vector_dist_celllist_hilb_performance_test )

///////////////////// INPUT DATA //////////////////////

// Cut-off radiuses. Can be put different number of values
openfpm::vector<float> r_cutoff {0.01, 0.02, 0.03};
// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The lower threshold for number of particles
size_t k_min = 15000;
// Ghost part of distributed vector
double ghost_part = 0.03;

///////////////////////////////////////////////////////

// Numbers of particles vector
openfpm::vector<size_t> n_particles;
// Vectors to store the data for 2D
openfpm::vector<openfpm::vector<double>> time_rand;
openfpm::vector<openfpm::vector<double>> time_hilb;
openfpm::vector<openfpm::vector<double>> time_total_rand;
openfpm::vector<openfpm::vector<double>> time_total_hilb;
// Vectors to store the data for 3D
openfpm::vector<openfpm::vector<double>> time_rand_2;
openfpm::vector<openfpm::vector<double>> time_hilb_2;
openfpm::vector<openfpm::vector<double>> time_total_rand_2;
openfpm::vector<openfpm::vector<double>> time_total_hilb_2;


BOOST_AUTO_TEST_CASE( vector_dist_celllist_random_test )
{
	//Benchmark test for 2D and 3D
	vd_celllist_random_benchmark<2>(k_start,k_min,ghost_part,r_cutoff,n_particles,time_rand,time_total_rand);
	vd_celllist_random_benchmark<3>(k_start,k_min,ghost_part,r_cutoff,n_particles,time_rand_2,time_total_rand_2);
}

BOOST_AUTO_TEST_CASE( vector_dist_celllist_hilbert_test )
{
	//Benchmark test for 2D and 3D
	vd_celllist_hilbert_benchmark<2>(k_start,k_min,ghost_part,r_cutoff,n_particles,time_hilb,time_total_hilb);
	vd_celllist_hilbert_benchmark<3>(k_start,k_min,ghost_part,r_cutoff,n_particles,time_hilb_2,time_total_hilb_2);
}

BOOST_AUTO_TEST_CASE(vector_dist_cl_performance_write_report)
{
	//Write report for 2D and 3D
	vd_celllist_performance_write_report<2>(r_cutoff,n_particles,time_hilb,time_rand,time_total_hilb,time_total_rand);
	vd_celllist_performance_write_report<3>(r_cutoff,n_particles,time_hilb_2,time_rand_2,time_total_hilb_2,time_total_rand_2);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_CL_HILB_PERFORMANCE_TESTS_HPP_ */
