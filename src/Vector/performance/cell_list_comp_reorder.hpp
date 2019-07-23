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
#include "cl_comp_performance_graph.hpp"

// Property tree
struct report_cell_list_func_tests
{
	boost::property_tree::ptree graphs;
};

report_cell_list_func_tests report_cl_funcs;

BOOST_AUTO_TEST_SUITE( celllist_getCellList_calc_forces_performance_test )


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

/*! \brief Function for random cell list test
 *
 */
template<unsigned int dim> void cell_list_getCellList_calc_force_benchmark(size_t cl_k_start,
		                                                                size_t cl_k_min,
																		openfpm::vector<float> & cl_r_cutoff,
																		openfpm::vector<size_t> & cl_n_particles,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_mean,
																		openfpm::vector<openfpm::vector<double>> & cl_time_rand_dev,
																		openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_mean,
																		openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_dev)
{
	report_cl_funcs.graphs.put("performance.celllist.dim",std::to_string(dim));

	cl_time_rand_mean.resize(cl_r_cutoff.size());
	cl_time_create_rand_mean.resize(cl_r_cutoff.size());
	cl_time_rand_dev.resize(cl_r_cutoff.size());
	cl_time_create_rand_dev.resize(cl_r_cutoff.size());

	std::string str("Testing " + std::to_string(dim) + "D vector, no order, cell-list");
	print_test_v(str,0);

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster<> & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			report_cl_funcs.graphs.put("performance.celllist.getCellList" + std::to_string(dim) + "D(" + std::to_string(r) + ").rcut",r_cut);

			int c = 0;

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with a random cell list k=" << k_int );

				report_cl_funcs.graphs.put("performance.celllist.getCellList" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").n",k_int);

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
				double sum_cl_mean = 0;
				double sum_cl_dev = 0;

				openfpm::vector<double> measures;
				for (size_t n = 0 ; n < N_STAT_TEST; n++)
				{
					measures.add(benchmark_get_celllist(NN,vd,r_cut));
				}
				standard_deviation(measures,sum_cl_mean,sum_cl_dev);
				//Average total time

				cl_time_create_rand_mean.get(r).add(sum_cl_mean);
				cl_time_create_rand_dev.get(r).add(sum_cl_dev);

				report_cl_funcs.graphs.put("performance.celllist.getCellList" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_cl_mean);
				report_cl_funcs.graphs.put("performance.celllist.getCellList" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_cl_dev);

				//Calculate forces

				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				measures.clear();
				for (size_t l = 0 ; l < N_STAT_TEST; l++)
				{measures.add(benchmark_calc_forces<dim>(NN,vd,r_cut));}
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				cl_time_rand_mean.get(r).add(sum_fr_mean);
				cl_time_rand_dev.get(r).add(sum_fr_dev);

				report_cl_funcs.graphs.put("performance.celllist.calc_forces" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_fr_mean);
				report_cl_funcs.graphs.put("performance.celllist.calc_forces" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_fr_dev);

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to create a cell-list: " << sum_cl_mean << " dev: " << sum_cl_dev << "    time to calculate forces: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;

				c++;
			}
		}
	}
}

/*! \brief Function for hilb cell list test
 *
 */
template<unsigned int dim> void cell_list_getCellList_hilb_calc_force_benchmark(size_t cl_k_start,
		                                                                 size_t cl_k_min,
																		 openfpm::vector<float> & cl_r_cutoff,
																		 openfpm::vector<size_t> & cl_n_particles,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_force_mean,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_force_dev,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_create_mean,
																		 openfpm::vector<openfpm::vector<double>> & cl_time_create_dev)
{
	report_cl_funcs.graphs.put("performance.celllist.dim",std::to_string(dim));

	cl_time_force_mean.resize(cl_r_cutoff.size());
	cl_time_create_mean.resize(cl_r_cutoff.size());
	cl_time_force_dev.resize(cl_r_cutoff.size());
	cl_time_create_dev.resize(cl_r_cutoff.size());

	std::string str("Testing " + std::to_string(dim) + "D vector, Hilbert comp reorder, cell list");
	print_test_v(str,0);

	{
		//For different r_cut
		for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
		{
			Vcluster<> & v_cl = create_vcluster();

			//Cut-off radius
			float r_cut = cl_r_cutoff.get(r);

			report_cl_funcs.graphs.put("performance.celllist.getCellList_hilb" + std::to_string(dim) + "D(" + std::to_string(r) + ").rcut",r_cut);

			//Number of particles
			size_t k = cl_k_start * v_cl.getProcessingUnits();

			//Counter number for amounts of particles
			size_t k_count = 1 + log2(k/cl_k_min);

			int c = 0;

			//For different number of particles
			for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
			{
				BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector with an Hilbert cell list k=" << k_int );

				report_cl_funcs.graphs.put("performance.celllist.getCellList_hilb" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").n",k_int);

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

				//Get a cell list hilb

				auto NN = vd.getCellList_hilb(r_cut);

				openfpm::vector<double> measures;

				double sum_cl_mean = 0;
				double sum_cl_dev = 0;
				for (size_t n = 0 ; n < N_STAT_TEST; n++)
				{measures.add(benchmark_get_celllist_hilb(NN,vd,r_cut));}
				standard_deviation(measures,sum_cl_mean,sum_cl_dev);
				//Average total time
				cl_time_create_mean.get(r).add(sum_cl_mean);
				cl_time_create_dev.get(r).add(sum_cl_dev);

				report_cl_funcs.graphs.put("performance.celllist.getCellList_hilb" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_cl_mean);
				report_cl_funcs.graphs.put("performance.celllist.getCellList_hilb" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_cl_dev);

				//Calculate forces

				double sum_fr_mean = 0;
				double sum_fr_dev = 0;

				// Initialize SFC (we are only interested in force calculation)
				NN.init_SFC();

				measures.clear();
				for (size_t l = 0 ; l < N_STAT_TEST; l++)
				{measures.add(benchmark_calc_forces_hilb<dim>(NN,vd,r_cut));}
				standard_deviation(measures,sum_fr_mean,sum_fr_dev);

				cl_time_force_mean.get(r).add(sum_fr_mean);
				cl_time_force_dev.get(r).add(sum_fr_dev);

				report_cl_funcs.graphs.put("performance.celllist.calc_forces_hilb" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").mean",sum_fr_mean);
				report_cl_funcs.graphs.put("performance.celllist.calc_forces_hilb" + std::to_string(dim) + "D(" + std::to_string(r) + ").npart(" + std::to_string(c) + ").dev",sum_fr_dev);


				if (v_cl.getProcessUnitID() == 0)
				{std::cout << "Cut-off = " << r_cut << ", Particles = " << k_int << ". Time to create: " << sum_cl_mean << " dev: " << sum_cl_dev << " time to calculate force: " << sum_fr_mean << " dev: " << sum_fr_dev << std::endl;}

				c++;
			}
		}
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_celllist_random_test )
{
	//Benchmark test for 2D and 3D
	cell_list_getCellList_calc_force_benchmark<3>(k_start,
			                                   k_min,
											   r_cutoff,
											   n_particles,
											   time_rand_mean,
											   time_rand_dev,
											   time_create_rand_mean,
											   time_create_rand_dev);


	cell_list_getCellList_calc_force_benchmark<2>(k_start,
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
	cell_list_getCellList_hilb_calc_force_benchmark<3>(k_start,
			                                    k_min,
												r_cutoff,
												n_particles,
												time_hilb_mean,
												time_hilb_dev,
												time_create_hilb_mean,
												time_create_hilb_dev);

	cell_list_getCellList_hilb_calc_force_benchmark<2>(k_start,
			                                    k_min,
												r_cutoff,
												n_particles,
												time_hilb_2_mean,
												time_hilb_2_dev,
												time_create_hilb_2_mean,
												time_create_hilb_2_dev);
}

/*! \brief Function for cell list hilb performance report
 *
 */
template<unsigned int dim> void cell_list_comp_reorder_report(GoogleChart & cg,
		                                                      openfpm::vector<float> & cl_r_cutoff,
															  openfpm::vector<size_t> & cl_n_particles,
															  openfpm::vector<openfpm::vector<double>> & cl_time_hilb_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_rand_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_hilb_dev,
															  openfpm::vector<openfpm::vector<double>> & cl_time_rand_dev,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_hilb_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_mean,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_hilb_dev,
															  openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_dev,
															  double & warning_level,
															  double & norm)
{
	cl_comp_normal_vs_hilbert_force_time<dim>(cg,
			                                  cl_n_particles,
											  cl_r_cutoff,
											  cl_time_hilb_mean,
											  cl_time_rand_mean,
											  cl_time_hilb_dev,
											  cl_time_rand_dev,
											  warning_level,
											  norm);

	cl_comp_normal_vs_hilbert_create_time<dim>(cg,
			                                   cl_n_particles,
											   cl_r_cutoff,
											   cl_time_create_hilb_mean,
											   cl_time_create_rand_mean,
											   cl_time_create_hilb_dev,
											   cl_time_create_rand_dev,
											   warning_level,
											   norm);
}

BOOST_AUTO_TEST_CASE(vector_dist_cl_performance_write_report)
{
	// Create a graphs

	//For different r_cut
/*	for (size_t r = 0; r < r_cutoff.size(); r++ )
	{
		report_cl_funcs.graphs.put("graphs.graph(" + std::to_string(r) + ").type","line");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").title","getCellList performance r_cut=" + std::to_string(r_cutoff.get(r)));
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").x.title","number of particles");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").y.title","Time seconds");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(0).source","performance.celllist.getCellList_hilb3D(" + std::to_string(r) + ").npart(#).mean");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").x.data(0).source","performance.celllist.getCellList_hilb3D(" + std::to_string(r) + ").npart(#).n");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(0).title","Cell-list hilbert");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(1).source","performance.celllist.getCellList3D(" + std::to_string(r) + ").npart(#).mean");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").x.data(1).source","performance.celllist.getCellList3D(" + std::to_string(r) + ").npart(#).n");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").y.data(1).title","Cell-list normal");
		report_cl_funcs.graphs.add("graphs.graph(" + std::to_string(r) + ").options.log_y","true");
	}

	boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
	boost::property_tree::write_xml("celllist_performance.xml", report_cl_funcs.graphs,std::locale(),settings);*/

	GoogleChart cg;

	std::string file_xml_ref(test_dir);
	file_xml_ref += std::string("/openfpm_pdata/celllist_performance_ref.xml");

	StandardXMLPerformanceGraph("celllist_performance.xml",file_xml_ref,cg);

	addUpdtateTime(cg);

	cg.write("celllist_performance.html");

/*
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
		cg.write("Celllist_comp_ord.html");
	}*/
}



BOOST_AUTO_TEST_SUITE_END()

#endif /* SRC_VECTOR_VECTOR_DIST_CL_HILB_PERFORMANCE_TESTS_HPP_ */
