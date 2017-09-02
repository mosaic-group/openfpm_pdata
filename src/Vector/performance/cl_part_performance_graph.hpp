/*
 * cl_part_performance_graph.hpp
 *
 *  Created on: Dec 23, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_PERFORMANCE_CL_PART_PERFORMANCE_GRAPH_HPP_
#define SRC_VECTOR_PERFORMANCE_CL_PART_PERFORMANCE_GRAPH_HPP_

#include "vector_dist_performance_util.hpp"

template<unsigned int dim> void cl_part_time(GoogleChart & cg,
		                                        openfpm::vector<size_t> & cl_n_particles,
												openfpm::vector<float> & cl_r_cutoff,
												openfpm::vector<size_t> &cl_orders,
												openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_hilb_mean,
												openfpm::vector<openfpm::vector<double>> cl_time_rand_mean,
												openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_hilb_dev,
												openfpm::vector<openfpm::vector<double>> cl_time_rand_dev)
{
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/cl_part_norm_hilbert_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/cl_part_norm_hilbert_dev_" + std::to_string(dim) + std::string("_ref"));

	std::string file_mean_save = std::string("cl_part_norm_hilbert_mean_" + std::to_string(dim) + std::to_string("_ref"));
	std::string file_var_save = std::string("cl_part_norm_hilbert_dev_" + std::to_string(dim) + std::to_string("_ref"));

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
			yp_mean.get(i).get(j).resize(1+cl_time_hilb_mean.get(i).get(j).size());
			yp_dev.get(i).get(j).resize(1+cl_time_hilb_dev.get(i).get(j).size());

			for (size_t k = 0 ; k < cl_time_hilb_mean.get(i).get(j).size() ; k++)
			{
				yp_mean.get(i).get(j).get(k) = cl_time_hilb_mean.get(i).get(j).get(k);
				yp_dev.get(i).get(j).get(k) = cl_time_hilb_dev.get(i).get(j).get(k);
			}
			yp_mean.get(i).get(j).get(cl_time_hilb_mean.get(i).get(j).size()) = cl_time_rand_mean.get(i).get(j);
			yp_dev.get(i).get(j).get(cl_time_hilb_mean.get(i).get(j).size()) = cl_time_rand_mean.get(i).get(j);
		}
	}

	names.add("No-order");
	for (size_t i = 0 ; i < cl_orders.size() ; i++)
		names.add(std::string("Order of: " + std::to_string(cl_orders.get(i))));

	for (size_t i = 0 ; i < cl_r_cutoff.size() ; i++)
		gnames.add("Cell-list performance, cut-off radius: " + std::to_string(cl_r_cutoff.get(i)));

	std::string y_string = std::string("Time to calculate forces (s)");
	std::string x_string = std::string("Number of particles");

	std::string str("<h1>Cell-list " + std::to_string(dim) + "-D performance tests: </h1>");
	str += "<h2> 1) Time to calculate forces in the case of particles randomly ordered in a vector of particles, and in the case of particles ordered along an hilbert curve of order N</h2>";

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
							 y_string,
							 false);
}

template<unsigned int dim> void cl_part_reorder_time(GoogleChart & cg,
		                                        openfpm::vector<size_t> & cl_n_particles,
												openfpm::vector<float> & cl_r_cutoff,
												openfpm::vector<size_t> &cl_orders,
												openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_reorder_mean,
												openfpm::vector<openfpm::vector<openfpm::vector<double>>> cl_time_reorder_dev)
{
	// Get the test dir
	std::string file_mean(test_dir);
	std::string file_var(test_dir);
	file_mean += std::string("/openfpm_pdata/cl_part_reorder_hilbert_mean_" + std::to_string(dim) + std::string("_ref"));
	file_var += std::string("/openfpm_pdata/cl_part_reorder_hilbert_dev_" + std::to_string(dim) + std::string("_ref"));

	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_mean;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_ref_dev;
	y_ref_mean.load(file_mean);
	y_ref_dev.load(file_var);

	// warning level
	openfpm::vector<int> warning_vlevel;

	// graphs data
	openfpm::vector<size_t> x;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y;
	openfpm::vector<openfpm::vector<openfpm::vector<double>>> y_dev;
	openfpm::vector<std::string> yn;

	for (size_t i = 0; i < cl_n_particles.size() ; i++)
		x.add(cl_n_particles.get(i));

	for (size_t i = 0; i < cl_orders.size(); i++)
		yn.add("Order of: " + std::to_string(cl_orders.get(i)));

	// Add reorder time
	y.resize(cl_time_reorder_mean.size());
	y_dev.resize(cl_time_reorder_dev.size());
	for (size_t r = 0; r < cl_time_reorder_mean.size(); r++)
	{
		y.get(r).resize(cl_time_reorder_mean.get(r).size());
		y_dev.get(r).resize(cl_time_reorder_mean.get(r).size());
		for (size_t k = 0; k < cl_time_reorder_mean.get(r).size(); k++)
		{
			// reorder time
			for (size_t m = 0; m < cl_time_reorder_mean.get(r).get(k).size(); m++)
			{
				// Put time
				y.get(r).get(k).add(cl_time_reorder_mean.get(r).get(k).get(m));
				y_dev.get(r).get(k).add(cl_time_reorder_dev.get(r).get(k).get(m));
			}
		}
	}

	// Save y
	y.save("cl_part_reorder_hilbert_mean_" + std::to_string(dim) + std::to_string("_ref"));
	y_dev.save("cl_part_reorder_hilbert_dev_" + std::to_string(dim) + std::to_string("_ref"));

	if (y_ref_mean.size() != 0)
	{
		yn.clear();
		y.clear();

		size_t i = cl_orders.size()-1;
		yn.add("Order of: " + std::to_string(cl_orders.get(i)));
		yn.add("interval");
		yn.add("interval");

		// Add reorder time
		y.resize(cl_time_reorder_mean.size());
		for (size_t r = 0; r < cl_time_reorder_mean.size(); r++)
		{
			int warning_level = -1;

			y.get(r).resize(cl_time_reorder_mean.get(r).size());
			for (size_t k = 0; k < cl_time_reorder_mean.get(r).size(); k++)
			{
				// reorder time
				size_t m = cl_orders.size()-1;

				// Put time
				y.get(r).get(k).add(cl_time_reorder_mean.get(r).get(k).get(m));
				y.get(r).get(k).add(y_ref_mean.get(r).get(k).get(m) - 3.0*y_ref_dev.get(r).get(k).get(m) );
				y.get(r).get(k).add(y_ref_mean.get(r).get(k).get(m) + 3.0*y_ref_dev.get(r).get(k).get(m) );

				warning_set(warning_level,cl_time_reorder_mean.get(r).get(k).get(m),y_ref_mean.get(r).get(k).get(m),y_ref_dev.get(r).get(k).get(m));
			}

			warning_vlevel.add(warning_level);
		}
	}

	// Speedup graphs report

	// Google charts options
	GCoptions options;

	options.yAxis = std::string("Time (s)");
	options.xAxis = std::string("Number of particles");
	options.lineWidth = 2;
	options.more = GC_ZOOM;

	//options.more = "hAxis: {logScale: true}";

	std::string str("<h1>Cell-list " + std::to_string(dim) + "-D performance tests: </h1>");
	str += "<h2> 1) Time to reorder the distributed vector</h2>";

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

#endif /* SRC_VECTOR_PERFORMANCE_CL_PART_PERFORMANCE_GRAPH_HPP_ */
