/*
 * vector_dist_performance_util.cpp
 *
 *  Created on: Feb 14, 2018
 *      Author: i-bird
 */

///////////////////////////// CONSTRUCT GRAPH //////////////////////////////

#include "vector_dist_performance_util.hpp"
#include "Plot/GoogleChart.hpp"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include "util/performance/performance_util.hpp"


/*! \brief Draw a standard performance graph
 *
 * \param file_mean
 *
 *
 */
void StandardPerformanceGraph(std::string file_mean,
		                      std::string file_var,
							  std::string file_mean_save,
							  std::string file_var_save,
							  GoogleChart & cg,
							  openfpm::vector<size_t> & xp,
							  openfpm::vector<openfpm::vector<openfpm::vector<double>>> & yp_mean,
							  openfpm::vector<openfpm::vector<openfpm::vector<double>>> & yp_dev,
							  openfpm::vector<std::string> & names,
							  openfpm::vector<std::string> & gnames,
							  std::string x_string,
							  std::string y_string,
							  bool use_log)
{
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

	if (names.size() == 0)
		return;

	for (size_t i = 0 ; i < names.size() ; i++)
		yn2.add(names.get(i));

	for (size_t i = 0; i < xp.size() ; i++)
		x.add(xp.get(i));

	yp_mean.save(file_mean_save);
	yp_dev.save(file_var_save);

	if (y_ref_mean.size() != 0 && yp_mean.size() != 0 && yp_mean.get(0).size() != 0)
	{
		// We reconstruct y and yn

		y2.clear();
		yn2.clear();

		for (size_t i = 0 ; i < yp_mean.get(0).get(0).size() ; i++)
		{
			yn2.add(names.get(i));
			yn2.add("interval");
			yn2.add("interval");
		}

		y2.resize(yp_mean.size());
		for (size_t r = 0; r < yp_mean.size(); r++)
		{
			int warning_level = -1;

			y2.get(r).resize(yp_mean.get(r).size());
			for (size_t k = 0; k < yp_mean.get(r).size(); k++)
			{

				// Number of graph points
				for (size_t g = 0 ; g < yp_mean.get(r).get(k).size() ; g++)
				{
					// Time for construction hilbert and random
					y2.get(r).get(k).add(yp_mean.get(r).get(k).get(g));
					y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(g) - 3.0*y_ref_dev.get(r).get(k).get(g));
					y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(g) + 3.0*y_ref_dev.get(r).get(k).get(g));

					warning_set(warning_level,yp_mean.get(r).get(k).get(g),y_ref_mean.get(r).get(k).get(g),y_ref_dev.get(r).get(k).get(g));
				}
			}

			warning_vlevel.add(warning_level);
		}
	}
	else
	{
		return;
	}

	// Calculation time graphs report

	// Google charts options
	GCoptions options2;

	options2.yAxis = std::string(y_string);
	options2.xAxis = std::string(x_string);
	options2.lineWidth = 4;

	for (size_t i = 0; i < y2.size() ; i++)
	{
		std::string chart_area;
		if (warning_vlevel.size() != 0)
			addchartarea(chart_area,warning_vlevel.get(i));

		if (use_log == true)
		{options2.more = GC_Y_LOG + "," + GC_ZOOM + chart_area;}
		else
		{options2.more = GC_ZOOM + chart_area;}

		options2.title = gnames.get(i);
		cg.AddLinesGraph(x,y2.get(i),yn2,options2);
	}
}


