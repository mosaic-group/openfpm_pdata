/*
 * performance_graph.hpp
 *
 *  Created on: Dec 23, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_PERFORMANCE_CL_COMP_PERFORMANCE_GRAPH_HPP_
#define SRC_VECTOR_PERFORMANCE_CL_COMP_PERFORMANCE_GRAPH_HPP_

#include "Plot/GoogleChart.hpp"

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
	// Get the test dir
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/cl_comp_norm_hilbert_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/cl_comp_norm_hilbert_dev_" + std::to_string(dim) + std::string("_ref"));

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_dev;
	y_ref_mean.load(file_mean);
	y_ref_dev.load(file_var);

	openfpm::vector<int> warning_vlevel;

	// time graphs data

	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_dev;

	openfpm::vector<std::string> yn;

	for (size_t i = 0; i < cl_n_particles.size() ; i++)
		x.add(cl_n_particles.get(i));

	yn.add("Random cell list");
	yn.add("Hilbert cell list");

	y.resize(cl_time_hilb_mean.size());
	y_dev.resize(cl_time_hilb_mean.size());
	for (size_t r = 0; r < cl_time_hilb_mean.size(); r++)
	{
		y.get(r).resize(cl_time_hilb_mean.get(r).size());
		y_dev.get(r).resize(cl_time_hilb_mean.get(r).size());
		for (size_t k = 0; k < cl_time_hilb_mean.get(r).size(); k++)
		{
			// Time for construction hilbert and random
			y.get(r).get(k).add(cl_time_rand_mean.get(r).get(k));
			y.get(r).get(k).add(cl_time_hilb_mean.get(r).get(k));

			y_dev.get(r).get(k).add(cl_time_rand_dev.get(r).get(k));
			y_dev.get(r).get(k).add(cl_time_hilb_dev.get(r).get(k));
		}
	}

	y.save("cl_comp_norm_hilbert_mean_" + std::to_string(dim) + std::to_string("_ref"));
	y_dev.save("cl_comp_norm_hilbert_dev_" + std::to_string(dim) + std::to_string("_ref"));

	if (y_ref_mean.size() != 0)
	{
		// We reconstruct y and yn

		y.clear();
		yn.clear();

		yn.add("Random cell list");
		yn.add("interval");
		yn.add("interval");
		yn.add("Hilbert cell list");
		yn.add("interval");
		yn.add("interval");
		y.resize(cl_time_hilb_mean.size());
		for (size_t r = 0; r < cl_time_hilb_mean.size(); r++)
		{
			int warning_level = -1;

			y.get(r).resize(cl_time_hilb_mean.get(r).size());
			for (size_t k = 0; k < cl_time_hilb_mean.get(r).size(); k++)
			{
				// Time for construction hilbert and random
				y.get(r).get(k).add(cl_time_rand_mean.get(r).get(k));

				y.get(r).get(k).add(y_ref_mean.get(r).get(k).get(0) - 3.0*y_ref_dev.get(r).get(k).get(0));
				y.get(r).get(k).add(y_ref_mean.get(r).get(k).get(0) + 3.0*y_ref_dev.get(r).get(k).get(0));
				y.get(r).get(k).add(cl_time_hilb_mean.get(r).get(k));
				y.get(r).get(k).add(y_ref_mean.get(r).get(k).get(1) - 3.0*y_ref_dev.get(r).get(k).get(1));
				y.get(r).get(k).add(y_ref_mean.get(r).get(k).get(1) + 3.0*y_ref_dev.get(r).get(k).get(1));

				warning_set(warning_level,cl_time_rand_mean.get(r).get(k),y_ref_mean.get(r).get(k).get(0),y_ref_dev.get(r).get(k).get(0));
				warning_set(warning_level,cl_time_hilb_mean.get(r).get(k),y_ref_mean.get(r).get(k).get(1),y_ref_dev.get(r).get(k).get(1));
			}

			warning_vlevel.add(warning_level);
		}
	}

	// Force time calculation

	// Google charts options
	GCoptions options;

	options.yAxis = std::string("Time to calculate force");
	options.xAxis = std::string("Number of particles");
	options.lineWidth = 4;

	std::string str("<h1>Cell-list " + std::to_string(dim) + "-D performance test: </h1>");
	str += "<h2> 1) Time to calculate forces</h2>";
	cg.addHTML(str);

	for (size_t i = 0; i < cl_r_cutoff.size(); i++)
	{
		std::string chart_area;
		if (warning_vlevel.size() != 0)
			addchartarea(chart_area,warning_vlevel.get(i));
		options.more = GC_Y_LOG + "," + GC_ZOOM + chart_area;

		options.title = std::string("Cell-list performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));
		cg.AddLinesGraph(x,y.get(i),yn,options);
	}
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
	// Get the test dir
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/cl_comp_create_norm_hilbert_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/cl_comp_create_norm_hilbert_dev_" + std::to_string(dim) + std::string("_ref"));

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_dev;
	y_ref_mean.load(file_mean);
	y_ref_dev.load(file_var);

	// warning level
	openfpm::vector<int> warning_vlevel;

	// Calculation time graphs data

	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y2;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y2_dev;
	openfpm::vector<std::string> yn2;

	yn2.add("Random cell list");
	yn2.add("Hilbert cell list");

	for (size_t i = 0; i < cl_n_particles.size() ; i++)
		x.add(cl_n_particles.get(i));

	y2.resize(cl_time_create_hilb_mean.size());
	y2_dev.resize(cl_time_create_hilb_mean.size());
	for (size_t r = 0; r < cl_time_create_hilb_mean.size(); r++)
	{
		y2.get(r).resize(cl_time_create_hilb_mean.get(r).size());
		y2_dev.get(r).resize(cl_time_create_hilb_mean.get(r).size());
		for (size_t k = 0; k < cl_time_create_hilb_mean.get(r).size(); k++)
		{
			// Put a total time
			y2.get(r).get(k).add(cl_time_create_rand_mean.get(r).get(k));
			y2.get(r).get(k).add(cl_time_create_hilb_mean.get(r).get(k));

			y2_dev.get(r).get(k).add(cl_time_create_rand_dev.get(r).get(k));
			y2_dev.get(r).get(k).add(cl_time_create_hilb_dev.get(r).get(k));
		}
	}

	y2.save("cl_comp_create_norm_hilbert_mean_" + std::to_string(dim) + std::to_string("_ref"));
	y2_dev.save("cl_comp_create_norm_hilbert_dev_" + std::to_string(dim) + std::to_string("_ref"));

	if (y_ref_mean.size() != 0)
	{
		// We reconstruct y and yn

		y2.clear();
		yn2.clear();

		yn2.add("Random cell list");
		yn2.add("interval");
		yn2.add("interval");
		yn2.add("Hilbert cell list");
		yn2.add("interval");
		yn2.add("interval");

		y2.resize(cl_time_create_hilb_mean.size());
		for (size_t r = 0; r < cl_time_create_hilb_mean.size(); r++)
		{
			int warning_level = -1;

			y2.get(r).resize(cl_time_create_hilb_mean.get(r).size());
			for (size_t k = 0; k < cl_time_create_hilb_mean.get(r).size(); k++)
			{
				// Time for construction hilbert and random
				y2.get(r).get(k).add(cl_time_create_rand_mean.get(r).get(k));
				y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(0) - 3.0*y_ref_dev.get(r).get(k).get(0));
				y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(0) + 3.0*y_ref_dev.get(r).get(k).get(0));
				y2.get(r).get(k).add(cl_time_create_hilb_mean.get(r).get(k));
				y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(1) - 3.0*y_ref_dev.get(r).get(k).get(1));
				y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(1) + 3.0*y_ref_dev.get(r).get(k).get(1));

				warning_set(warning_level,cl_time_create_rand_mean.get(r).get(k),y_ref_mean.get(r).get(k).get(0),y_ref_dev.get(r).get(k).get(0));
				warning_set(warning_level,cl_time_create_rand_mean.get(r).get(k),y_ref_mean.get(r).get(k).get(1),y_ref_dev.get(r).get(k).get(1));
			}

			warning_vlevel.add(warning_level);
		}
	}

	// Calculation time graphs report

	// Google charts options
	GCoptions options2;

	options2.yAxis = std::string("Time to create the cell-list");
	options2.xAxis = std::string("Number of particles");
	options2.lineWidth = 4;

	std::string str2("<h2>2) Time to create the cell-list</h2>");
	cg.addHTML(str2);

	for (size_t i = 0; i < cl_r_cutoff.size(); i++)
	{
		std::string chart_area;
		if (warning_vlevel.size() != 0)
			addchartarea(chart_area,warning_vlevel.get(i));
		options2.more = GC_Y_LOG + "," + GC_ZOOM + chart_area;

		options2.title = std::string("Cell-list performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));
		cg.AddLinesGraph(x,y2.get(i),yn2,options2);
	}
}


///////////////////////////////// PARTICLE REORDERING //////////////////////////////////////////////////////



#endif /* SRC_VECTOR_PERFORMANCE_CL_COMP_PERFORMANCE_GRAPH_HPP_ */
