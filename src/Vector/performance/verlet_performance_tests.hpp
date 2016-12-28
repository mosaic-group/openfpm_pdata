/*
 * vector_dist_verlet_performance_tests.hpp
 *
 *  Created on: Mar 9, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_VERLET_PERFORMANCE_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_VERLET_PERFORMANCE_TESTS_HPP_

BOOST_AUTO_TEST_SUITE( verletlist_part_reorder_performance_test )

///////////////////// INPUT DATA //////////////////////

// Cut-off radiuses. Can be put different number of values
openfpm::vector<float> r_cutoff {0.004, 0.007, 0.01};
// Orders of a curve. Can be put different number of values
// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The minimal amount of particles
size_t k_min = 15000;
// Ghost part of distributed vector

///////////////////////////////////////////////////////

// Numbers of particles vector
openfpm::vector<size_t> n_particles;
// Vectors to store the data for 2D
openfpm::vector<openfpm::vector<double>> time_force_mean;
openfpm::vector<openfpm::vector<double>> time_force_dev;
openfpm::vector<openfpm::vector<double>> time_create_mean;
openfpm::vector<openfpm::vector<double>> time_create_dev;

// Vectors to store the data for 3D
openfpm::vector<openfpm::vector<double>> time_force_mean_2;
openfpm::vector<openfpm::vector<double>> time_force_dev_2;
openfpm::vector<openfpm::vector<double>> time_create_mean_2;
openfpm::vector<openfpm::vector<double>> time_create_dev_2;


BOOST_AUTO_TEST_CASE( vector_dist_verlet_test )
{
	//Benchmark test for 2D and 3D
	vd_verlet_random_benchmark<3>(k_start,k_min,r_cutoff,n_particles,time_force_mean,time_create_mean,time_force_dev,time_create_dev);
	vd_verlet_random_benchmark<2>(k_start,k_min,r_cutoff,n_particles,time_force_mean_2,time_create_mean_2,time_force_dev_2,time_create_dev_2);
}

BOOST_AUTO_TEST_CASE(vector_dist_verlet_performance_write_report)
{
	GoogleChart cg;

	//Write report for 2D and 3D
	vd_verlet_performance_write_report<3>(cg,r_cutoff,n_particles,time_force_mean,time_force_dev,time_create_mean,time_create_dev);
	vd_verlet_performance_write_report<2>(cg,r_cutoff,n_particles,time_force_mean_2,time_force_dev_2,time_create_mean_2,time_create_dev_2);

	if (create_vcluster().getProcessUnitID() == 0)
	{
		addUpdtateTime(cg);

		cg.write("Verletlist_comp.html");
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_VERLET_PERFORMANCE_TESTS_HPP_ */
