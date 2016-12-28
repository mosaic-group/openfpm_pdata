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

BOOST_AUTO_TEST_SUITE( celllist_comp_reorder_performance_test )

///////////////////// INPUT DATA //////////////////////

// Cut-off radiuses. Can be put different number of values
openfpm::vector<float> r_cutoff {0.004, 0.007, 0.01};
// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
size_t k_start = 100000;
// The lower threshold for number of particles
size_t k_min = 15000;

///////////////////////////////////////////////////////

// Numbers of particles vector
openfpm::vector<size_t> n_particles;
// Vectors to store the data for 2D
openfpm::vector<openfpm::vector<double>> time_rand_mean;
openfpm::vector<openfpm::vector<double>> time_hilb_mean;
openfpm::vector<openfpm::vector<double>> time_rand_dev;
openfpm::vector<openfpm::vector<double>> time_hilb_dev;
openfpm::vector<openfpm::vector<double>> time_create_rand_mean;
openfpm::vector<openfpm::vector<double>> time_create_hilb_mean;
openfpm::vector<openfpm::vector<double>> time_create_rand_dev;
openfpm::vector<openfpm::vector<double>> time_create_hilb_dev;
// Vectors to store the data for 3D
openfpm::vector<openfpm::vector<double>> time_rand_2_mean;
openfpm::vector<openfpm::vector<double>> time_hilb_2_mean;
openfpm::vector<openfpm::vector<double>> time_rand_2_dev;
openfpm::vector<openfpm::vector<double>> time_hilb_2_dev;
openfpm::vector<openfpm::vector<double>> time_create_rand_2_mean;
openfpm::vector<openfpm::vector<double>> time_create_hilb_2_mean;
openfpm::vector<openfpm::vector<double>> time_create_rand_2_dev;
openfpm::vector<openfpm::vector<double>> time_create_hilb_2_dev;


BOOST_AUTO_TEST_CASE( vector_dist_celllist_random_test )
{
	//Benchmark test for 2D and 3D
	cell_list_comp_reorder_random_benchmark<3>(k_start,
			                                   k_min,
											   r_cutoff,
											   n_particles,
											   time_rand_mean,
											   time_rand_dev,
											   time_create_rand_mean,
											   time_create_rand_dev);


	cell_list_comp_reorder_random_benchmark<2>(k_start,
			                                   k_min,
											   r_cutoff,
											   n_particles,
											   time_rand_2_mean,
											   time_rand_2_dev,
											   time_create_rand_2_mean,
											   time_create_rand_2_dev);
}

BOOST_AUTO_TEST_CASE( vector_dist_celllist_hilbert_test )
{
	//Benchmark test for 2D and 3D
	cell_list_comp_reorder_hilbert_benchmark<3>(k_start,
			                                    k_min,
												r_cutoff,
												n_particles,
												time_hilb_mean,
												time_hilb_dev,
												time_create_hilb_mean,
												time_create_hilb_dev);

	cell_list_comp_reorder_hilbert_benchmark<2>(k_start,
			                                    k_min,
												r_cutoff,
												n_particles,
												time_hilb_2_mean,
												time_hilb_2_dev,
												time_create_hilb_2_mean,
												time_create_hilb_2_dev);
}

BOOST_AUTO_TEST_CASE(vector_dist_cl_performance_write_report)
{
	GoogleChart cg;

	double warning_level = 0;
	double norm = 0;

	//Write report for 2D and 3D
	cell_list_comp_reorder_report<3>(cg,
			                         r_cutoff,
									 n_particles,
									 time_hilb_mean,
									 time_rand_mean,
									 time_hilb_dev,
									 time_rand_dev,
									 time_create_hilb_mean,
									 time_create_rand_mean,
									 time_create_hilb_dev,
									 time_create_rand_dev,
									 warning_level,
									 norm);

	cell_list_comp_reorder_report<2>(cg,
			                         r_cutoff,
									 n_particles,
									 time_hilb_2_mean,
									 time_rand_2_mean,
									 time_rand_2_dev,
									 time_hilb_2_dev,
									 time_create_hilb_2_mean,
									 time_create_rand_2_mean,
									 time_create_hilb_2_dev,
									 time_create_rand_2_dev,
									 warning_level,
									 norm);

	addUpdtateTime(cg);

	if (create_vcluster().getProcessUnitID() == 0)
	{
		// write the xml report
		pt.put("celllist.comp.warning",warning_level);

		cg.write(std::string(test_dir) + "/openfpm_pdata/Celllist_comp_ord.html");
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_CL_HILB_PERFORMANCE_TESTS_HPP_ */
