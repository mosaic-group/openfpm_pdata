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
#include "Plot/GoogleChart.hpp"
#include <functional>

BOOST_AUTO_TEST_SUITE( celllist_part_reorder_performance_test )

///////////////////// INPUT DATA //////////////////////

// Cut-off radiuses. Can be put different number of values
openfpm::vector<float> r_cutoff {0.004, 0.007, 0.01};
// Orders of a curve. Can be put different number of values
openfpm::vector<size_t> orders = {1,2,3};
// Number of steps of moving the particles
size_t n_moving = 8;
// Moving distance (step size)
double dist = 0.03;
// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The minimal amount of particles
size_t k_min = 15000;

///////////////////////////////////////////////////////

// Numbers of particles vector
openfpm::vector<size_t> n_particles;
// Vectors to store the data for 2D
openfpm::vector<openfpm::vector<double>> time_rand_mean;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_hilb_mean;
openfpm::vector<openfpm::vector<double>> time_rand_dev;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_hilb_dev;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_reorder_mean;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_reorder_dev;
// Vectors to store the data for 3D
openfpm::vector<openfpm::vector<double>> time_rand_2_mean;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_hilb_2_mean;
openfpm::vector<openfpm::vector<double>> time_rand_2_dev;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_hilb_2_dev;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_reorder_2_mean;
openfpm::vector<openfpm::vector<openfpm::vector<double>>> time_reorder_2_dev;


BOOST_AUTO_TEST_CASE( vector_dist_cl_random_test )
{
	//Benchmark test for 2D and 3D
	cell_list_part_reorder_random_benchmark<3>(k_start,k_min,r_cutoff,n_particles,time_rand_mean,time_rand_dev);
	cell_list_part_reorder_random_benchmark<2>(k_start,k_min,r_cutoff,n_particles,time_rand_2_mean,time_rand_2_dev);
}

BOOST_AUTO_TEST_CASE( vector_dist_cl_hilbert_test )
{
	//Benchmark test for 2D and 3D
	cell_list_part_reorder_hilbert_benchmark<3>(k_start,
			                                    k_min,
												n_moving,
												dist,
												r_cutoff,
												n_particles,
												orders,
												time_hilb_mean,
												time_reorder_mean,
												time_hilb_dev,
												time_reorder_dev);

	cell_list_part_reorder_hilbert_benchmark<2>(k_start,
			                                    k_min,
												n_moving,
												dist,
												r_cutoff,
												n_particles,
												orders,
												time_hilb_2_mean,
												time_reorder_2_mean,
												time_hilb_2_dev,
												time_reorder_2_dev);
}

BOOST_AUTO_TEST_CASE(vector_dist_cl_performance_write_report)
{
	GoogleChart cg;

/*	n_particles.save("debug1");
	time_rand_mean.save("debug2");
	time_hilb_mean.save("debug3");
	time_rand_dev.save("debug4");
	time_hilb_dev.save("debug5");
	time_reorder_mean.save("debug6");
	time_reorder_dev.save("debug7");
	time_hilb_moved.save("debug8");

	time_rand_2_mean.save("debug9");
	time_hilb_2_mean.save("debug10");
	time_rand_2_dev.save("debug11");
	time_hilb_2_dev.save("debug12");
	time_reorder_2_mean.save("debug13");
	time_reorder_2_dev.save("debug14");
	time_hilb_moved_2.save("debug15");*/

/*	n_particles.load("debug1");
	time_rand_mean.load("debug2");
	time_hilb_mean.load("debug3");
	time_rand_dev.load("debug4");
	time_hilb_dev.load("debug5");
	time_reorder_mean.load("debug6");
	time_reorder_dev.load("debug7");
	time_hilb_moved.load("debug8");

	time_rand_2_mean.load("debug9");
	time_hilb_2_mean.load("debug10");
	time_rand_2_dev.load("debug11");
	time_hilb_2_dev.load("debug12");
	time_reorder_2_mean.load("debug13");
	time_reorder_2_dev.load("debug14");
	time_hilb_moved_2.load("debug15");*/

	//Write report for 2D and 3D
	cell_list_part_reorder_report<3>(cg,
			                         n_moving,
			                         r_cutoff,
									 n_particles,
									 orders,
									 time_hilb_mean,
									 time_rand_mean,
									 time_reorder_mean,
									 time_hilb_dev,
									 time_rand_dev,
									 time_reorder_dev);

	cell_list_part_reorder_report<2>(cg,
			                         n_moving,
			                         r_cutoff,
									 n_particles,
									 orders,
									 time_hilb_2_mean,
									 time_rand_2_mean,
									 time_reorder_2_mean,
									 time_hilb_2_dev,
									 time_rand_2_dev,
									 time_reorder_2_dev);

	if (create_vcluster().getProcessUnitID() == 0)
		cg.write(std::string(test_dir) + "/openfpm_pdata/Celllist_part_ord.html");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_ */
