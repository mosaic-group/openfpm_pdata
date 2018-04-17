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
#include "cl_part_performance_graph.hpp"
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


/*! \brief Function for cell list test without an Hilbert curve reordering (unordered positioning)
 *
 */
template<unsigned int dim> void cell_list_part_reorder_random_benchmark(size_t cl_k_start,
		                                                                size_t cl_k_min,
																		openfpm::vector<float> & cl_r_cutoff,
																		openfpm::vector<size_t> & cl_n_particles,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_mean,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_dev)
{
	cl_time_rand_mean.resize(cl_r_cutoff.size());
	cl_time_rand_dev.resize(cl_r_cutoff.size());

	std::string str("Testing " + std::to_string(dim) + "D vector, no-order, Cell-list");
	print_test_v(str,0);

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector without an Hilbert curve reordering k=" << k_int );

				if (cl_n_particles.size() < k_count)
					cl_n_particles.add(k_int);

				Box<dim,float> box;

				for (size_t i = 0; i < dim; i++)
				{
					box.setLow(i,0.0);
					box.setHigh(i,1.0);
				}

				// Boundary conditions
				size_t bc[dim];

				for (size_t i = 0; i < dim; i++)
					bc[i] = PERIODIC;

				vector_dist<dim,float, aggregate<float[dim]> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

				// Initialize a dist vector
				vd_initialize<dim>(vd, v_cl, k_int);

				vd.template ghost_get<0>();

				//Get a cell list

				auto NN = vd.getCellList(r_cut);
				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				benchmark_get_celllist(NN,vd,r_cut);

				//Calculate forces
				size_t l = 0;

				openfpm::vector<double> measures;
				for ( ; l < N_STAT_TEST; l++)
					measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				cl_time_rand_mean.get(r).add(sum_fr_mean);
				cl_time_rand_dev.get(r).add(sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << " time to calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;
			}
		}
	}
}


/*! \brief Function for cell list test with an Hilbert curve reordering
 *
 */
template<unsigned int dim> void cell_list_part_reorder_hilbert_benchmark(size_t cl_k_start,
		                                                                 size_t cl_k_min,
																		 size_t n_moving,
																		 double dist,
																		 openfpm::vector<float> & cl_r_cutoff,
																		 openfpm::vector<size_t> & cl_n_particles,
																		 openfpm::vector<size_t> &cl_orders,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_hilb_mean,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_reorder_hilb_mean,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_hilb_dev,
																		 openfpm::vector<openfpm::vector<openfpm::vector<double>>> &cl_time_reorder_hilb_dev)
{
	{
		cl_time_hilb_mean.resize(cl_r_cutoff.size());
		cl_time_hilb_dev.resize(cl_r_cutoff.size());
		for (size_t r = 0; r < cl_time_hilb_mean.size(); r++)
		{
			cl_time_hilb_mean.get(r).resize(cl_n_particles.size());
			cl_time_hilb_dev.get(r).resize(cl_n_particles.size());
			for (size_t k = 0; k < cl_time_hilb_mean.get(r).size(); k++)
			{
				cl_time_hilb_mean.get(r).get(k).resize(cl_orders.size());
				cl_time_hilb_dev.get(r).get(k).resize(cl_orders.size());
			}
		}


		cl_time_reorder_hilb_mean.resize(cl_r_cutoff.size());
		cl_time_reorder_hilb_dev.resize(cl_r_cutoff.size());
		for (size_t r = 0; r < cl_time_reorder_hilb_mean.size(); r++)
		{
			cl_time_reorder_hilb_mean.get(r).resize(cl_n_particles.size());
			cl_time_reorder_hilb_dev.get(r).resize(cl_n_particles.size());
			for (size_t k = 0; k < cl_time_reorder_hilb_mean.get(r).size(); k++)
			{
				cl_time_reorder_hilb_mean.get(r).get(k).resize(cl_orders.size());
				cl_time_reorder_hilb_dev.get(r).get(k).resize(cl_orders.size());
			}
		}

		// Print test
		std::string str("Testing " + std::to_string(dim) + "D vector, Hilbert curve reordering, Cell-List");
		print_test_v(str,0);

		// For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster & v_cl = create_vcluster();

			// Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			// Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//For different curve orders
			for ( size_t i = 0; i < cl_orders.size(); i++)
			{
				size_t m = cl_orders.get(i);
				size_t part = 0;

				for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2, part++ )
				{
					BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert curve reordering k=" << k_int );

					Box<dim,float> box;

					for (size_t i = 0; i < dim; i++)
					{
						box.setLow(i,0.0);
						box.setHigh(i,1.0);
					}

					// Boundary conditions
					size_t bc[dim];

					for (size_t i = 0; i < dim; i++)
						bc[i] = PERIODIC;

					vector_dist<dim,float, aggregate<float[dim]> > vd(k_int,box,bc,Ghost<dim,float>(r_cut));

					// Initialize a dist vector
					vd_initialize<dim>(vd, v_cl, k_int);

					//Reorder a vector

					double sum_reorder_mean = 0;
					double sum_reorder_dev = 0;

					openfpm::vector<double> measures;

					for (size_t h = 0 ; h < N_STAT_TEST; h++)
						measures.add(benchmark_reorder(vd,m));
					standard_deviation(measures,sum_reorder_mean,sum_reorder_dev);

					//Average reorder time
					cl_time_reorder_hilb_mean.get(r).get(part).get(i) = sum_reorder_mean;
					cl_time_reorder_hilb_dev.get(r).get(part).get(i) = sum_reorder_dev;

					vd.template ghost_get<0>();

					//Get cell list

					auto NN = vd.getCellList(r_cut);
					benchmark_get_celllist(NN,vd,r_cut);

					//Calculate forces

					double sum_fr_mean = 0;
					double sum_fr_dev = 0;
					measures.clear();

					for (size_t l = 0 ; l < N_STAT_TEST ; l++)
						measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));
					standard_deviation(measures,sum_fr_mean,sum_fr_dev);

					cl_time_hilb_mean.get(r).get(part).get(i) = sum_fr_mean;
					cl_time_hilb_dev.get(r).get(part).get(i) = sum_fr_dev;


					if (v_cl.getProcessUnitID() == 0)
						std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to reorder: " << sum_reorder_mean << " dev: " << sum_reorder_dev << "      time calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;
				}
			}
		}
	}
}



/*! \brief Function for cell list performance report
 *
 */
template<unsigned int dim> void cell_list_part_reorder_report(GoogleChart & cg,
		                                                      size_t n_moving,
		                                                      openfpm::vector<float> & cl_r_cutoff,
															  openfpm::vector<size_t> & cl_n_particles,
															  openfpm::vector<size_t> cl_orders,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_hilb_mean,
															  openfpm::vector<openfpm::vector<double>> cl_time_rand_mean,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_reorder_mean,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_hilb_dev,
															  openfpm::vector<openfpm::vector<double>> cl_time_rand_dev,
															  openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_reorder_dev)
{
	openfpm::vector<size_t> x;

	for (size_t i = 0; i < cl_n_particles.size() ; i++)
		x.add(cl_n_particles.get(i));

	// Speedup graphs data

	cl_part_time<dim>(cg,
			          cl_n_particles,
					  cl_r_cutoff,
					  cl_orders,
					  cl_time_hilb_mean,
					  cl_time_rand_mean,
					  cl_time_hilb_dev,
					  cl_time_rand_dev);

	cl_part_reorder_time<dim>(cg,
			                  cl_n_particles,
							  cl_r_cutoff,
							  cl_orders,
							  cl_time_reorder_mean,
							  cl_time_reorder_dev);
}



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

	addUpdtateTime(cg);

	if (create_vcluster().getProcessUnitID() == 0)
		cg.write("Celllist_part_ord.html");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_CL_PERFORMANCE_TESTS_HPP_ */
