/*
 * performance_graph.hpp
 *
 *  Created on: Dec 23, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_PERFORMANCE_CL_COMP_PERFORMANCE_GRAPH_HPP_
#define SRC_VECTOR_PERFORMANCE_CL_COMP_PERFORMANCE_GRAPH_HPP_

#include "Plot/GoogleChart.hpp"
#include "vector_dist_performance_util.hpp"

/////////////////////////// COMPUTATIONAL HILBERT CURVE ORDERING ///////////////////////////////////////

/*! \brief Time to calculate forces using random order or hilber curve compute ordering
 *
 * \tparam dim dimensionality of the test
 *
 * \param cg GoogleChart object
 * \param cl_n_particles number of particles for each test
 * \param cl_r_cutoff cell-list spacing used to construct a cell-list
 * \param cl_time_hilb for each particle set, for each cut-off radius time to calculate the forces (Hilbert curve order)
 * \param cl_time_rand for each particle set, for each cut-off radius, time to calculate the forces (no order)
 *
 */
template<unsigned int dim> void cl_comp_normal_vs_hilbert_force_time(GoogleChart & cg,
		                                                             openfpm::vector<size_t> & cl_n_particles,
																	 openfpm::vector<float> & cl_r_cutoff,
																	 openfpm::vector<openfpm::vector<double>> & cl_time_hilb_mean,
																	 openfpm::vector<openfpm::vector<double>> & cl_time_rand_mean,
																	 openfpm::vector<openfpm::vector<double>> & cl_time_hilb_dev,
																	 openfpm::vector<openfpm::vector<double>> & cl_time_rand_dev,
																	 double & warning_level,
																	 double & norm)
{
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/cl_comp_norm_hilbert_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/cl_comp_norm_hilbert_dev_" + std::to_string(dim) + std::string("_ref"));

	std::string file_mean_save = std::string("cl_comp_norm_hilbert_mean_" + std::to_string(dim) + std::to_string("_ref"));
	std::string file_var_save = std::string("cl_comp_norm_hilbert_dev_" + std::to_string(dim) + std::to_string("_ref"));

	openfpm::vector<size_t> xp = cl_n_particles;

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_dev;

	openfpm::vector<std::string> names;
	openfpm::vector<std::string> gnames;

	yp_mean.resize(cl_time_rand_mean.size());
	yp_dev.resize(cl_time_rand_dev.size());
	for (size_t i = 0 ; i < yp_mean.size() ; i++)
	{
		yp_mean.get(i).resize(cl_time_rand_mean.get(i).size());
		yp_dev.get(i).resize(cl_time_rand_dev.get(i).size());

		for (size_t j = 0 ; j < yp_mean.get(i).size() ; j++)
		{
			yp_mean.get(i).get(j).resize(2);
			yp_dev.get(i).get(j).resize(2);

			yp_mean.get(i).get(j).get(0) = cl_time_hilb_mean.get(i).get(j);
			yp_dev.get(i).get(j).get(0) = cl_time_hilb_dev.get(i).get(j);

			yp_mean.get(i).get(j).get(1) = cl_time_rand_mean.get(i).get(j);
			yp_dev.get(i).get(j).get(1) = cl_time_rand_mean.get(i).get(j);
		}
	}

	names.add("Random cell list");
	names.add("Hilbert cell list");

	for (size_t i = 0 ; i < cl_r_cutoff.size() ; i++)
		gnames.add("Cell-list performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));

	std::string y_string = std::string("Time to calculate forces");
	std::string x_string = std::string("Number of particles");

	std::string str("<h1>Cell-list " + std::to_string(dim) + "-D performance test: </h1>");
	str += "<h2> 1) Time to calculate forces</h2>";
	cg.addHTML(str);

	StandardPerformanceGraph(file_mean,
			                 file_var,
							 file_mean_save,
							 file_var_save,
							 cg,
							 xp,
							 yp_mean,
							 yp_dev,
							 names,
							 gnames,
							 x_string,
							 y_string);
}

/*! \brief Output the graph normal cell-list vs Hilbert cell-list (Total time)
 *
 * \tparam dim dimensionality of the test
 *
 * \param cg GoogleChart object
 * \param cl_n_particles number of particles for each test
 * \param cl_r_cutoff cell-list spacing used to construct a cell-list
 * \param cl_time_force_hilb for each particle set, for each cut-off radius, for each order time to calculate the forces
 * \param cl_time_force_rand for each particle set, for each cut-off radius, time to calculate the forces
 *
 */
template<unsigned int dim> void cl_comp_normal_vs_hilbert_create_time(GoogleChart & cg,
		                                                              openfpm::vector<size_t> & cl_n_particles,
																	  openfpm::vector<float> & cl_r_cutoff,
																	  openfpm::vector<openfpm::vector<double>> & cl_time_create_hilb_mean,
																	  openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_mean,
																	  openfpm::vector<openfpm::vector<double>> & cl_time_create_hilb_dev,
																	  openfpm::vector<openfpm::vector<double>> & cl_time_create_rand_dev,
																	  double & warning_level,
																	  double & norm)
{
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/cl_comp_create_norm_hilbert_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/cl_comp_create_norm_hilbert_dev_" + std::to_string(dim) + std::string("_ref"));

	std::string file_mean_save = std::string("cl_comp_create_norm_hilbert_mean_" + std::to_string(dim) + std::string("_ref"));
	std::string file_var_save = std::string("cl_comp_create_norm_hilbert_dev_" + std::to_string(dim) + std::string("_ref"));

	openfpm::vector<size_t> xp = cl_n_particles;

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> yp_dev;

	openfpm::vector<std::string> names;
	openfpm::vector<std::string> gnames;

	yp_mean.resize(cl_time_create_rand_mean.size());
	yp_dev.resize(cl_time_create_rand_dev.size());
	for (size_t i = 0 ; i < yp_mean.size() ; i++)
	{
		yp_mean.get(i).resize(cl_time_create_rand_mean.get(i).size());
		yp_dev.get(i).resize(cl_time_create_rand_dev.get(i).size());

		for (size_t j = 0 ; j < yp_mean.get(i).size() ; j++)
		{
			yp_mean.get(i).get(j).resize(2);
			yp_dev.get(i).get(j).resize(2);

			yp_mean.get(i).get(j).get(0) = cl_time_create_hilb_mean.get(i).get(j);
			yp_dev.get(i).get(j).get(0) = cl_time_create_hilb_dev.get(i).get(j);

			yp_mean.get(i).get(j).get(1) = cl_time_create_rand_mean.get(i).get(j);
			yp_dev.get(i).get(j).get(1) = cl_time_create_rand_mean.get(i).get(j);
		}
	}

	names.add("Random cell list");
	names.add("Hilbert cell list");

	for (size_t i = 0 ; i < cl_r_cutoff.size() ; i++)
		gnames.add("Cell-list performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));

	std::string y_string = std::string("Time to create the cell-list");
	std::string x_string = std::string("Number of particles");

	std::string str("<h1>Cell-list " + std::to_string(dim) + "-D performance test: </h1>");
	str += "<h2> 1) Time to create the cell-list</h2>";
	cg.addHTML(str);

	StandardPerformanceGraph(file_mean,
			                 file_var,
							 file_mean_save,
							 file_var_save,
							 cg,
							 xp,
							 yp_mean,
							 yp_dev,
							 names,
							 gnames,
							 x_string,
							 y_string);
}


///////////////////////////////// PARTICLE REORDERING //////////////////////////////////////////////////////



#endif /* SRC_VECTOR_PERFORMANCE_CL_COMP_PERFORMANCE_GRAPH_HPP_ */
