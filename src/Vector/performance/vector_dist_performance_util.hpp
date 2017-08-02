/*
 * vector_dist_performance_util.hpp
 *
 *  Created on: Jun 17, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_
#define SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_

#include "Plot/GoogleChart.hpp"
#include "vector_dist_performance_common.hpp"

// Number of tests
#define N_VERLET_TEST 3
#define N_STAT_TEST 30

static inline void warning_set(int & warning_level, double mean, double mean_ref, double sigma)
{
	int warning_level_candidate;

	if (mean - mean_ref < -2.0*sigma )
		warning_level_candidate = -1;
	else if (mean - mean_ref < 2.0*sigma)
		warning_level_candidate = 0;
	else if (mean - mean_ref < 3.0*sigma)
		warning_level_candidate = 1;
	else
		warning_level_candidate = 2;

	if (warning_level_candidate > warning_level)
		warning_level = warning_level_candidate;
}

static inline void addchartarea(std::string & chart_area, int lvl)
{
	std::string color;

	if (lvl == -1)
	{
		chart_area = std::string(",chartArea: {\
		    backgroundColor: {\
		        stroke: '#00FF00',\
		        strokeWidth: 6\
		    }\
		}");
	}
	else if (lvl == 0)
	{
		// NOTHING TO DO
	}
	else if (lvl == 1)
	{
		chart_area = std::string(",chartArea: {\
		    backgroundColor: {\
		        stroke: '#FFFF00',\
		        strokeWidth: 6\
		    }\
		}");
	}
	else if (lvl == 2)
	{
		chart_area = std::string(",chartArea: {\
		    backgroundColor: {\
		        stroke: '#FF0000',\
		        strokeWidth: 6\
		    }\
		}");
	}

}

void addUpdtateTime(GoogleChart & cg)
{
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );

    std::stringstream str;

    str << "<h3>Updated: " << now->tm_mday << "/" << now->tm_mon + 1 << "/" << now->tm_year+1900 << "     " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << std::endl;

	cg.addHTML(str.str());
}

/*! \brief Standard deviation
 *
 * \param measures set of measures
 * \param mean the mean of the measures
 *
 * \return the standard deviation
 *
 */
static inline void standard_deviation(openfpm::vector<double> measures, double & mean, double & dev)
{
	for (size_t i = 0 ; i < measures.size() ; i++)
		mean += measures.get(i);
	mean /= measures.size();

	dev = 0;
	for (size_t i = 0 ; i < measures.size() ; i++)
		dev += (measures.get(i) - mean)*(measures.get(i) - mean);

	dev = sqrt(dev / (measures.size() - 1));
}

/*! \brief Print out only ones (no matter how many processors involved)
 *
 * \param test, sz Data to print out
 */
void print_test_v(std::string test)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test  << "\n";
}


/*! \brief Benchmark particles' forces time
 *
 * \param NN Cell list
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> double benchmark_calc_forces(T & NN, V & vd, float r_cut)
{
	//Timer
	timer t;
	t.start();

	calc_forces<dim,prp>(NN,vd,r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark reordering time
 *
 * \param vd Distributed vector
 * \param m Order of an Hilbert curve
 *
 * \return real time
 */
template<typename V> double benchmark_reorder(V & vd, size_t m)
{
	//Timer
	timer t;
	t.start();

	//Reorder
	vd.reorder(m);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark celllist getting time
 *
 * \param NN Cell list
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<typename T, typename V> double benchmark_get_celllist(T & NN, V & vd, float r_cut)
{
	// Cell list timer
	timer t;
	t.start();

	//get cell list
	NN = vd.getCellList(r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark verlet getting time
 *
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<typename V> double benchmark_get_verlet(V & vd, float r_cut)
{
	openfpm::vector<openfpm::vector<size_t>> verlet;
	//Timer
	timer t;
	t.start();

	//get verlet
	auto vr = vd.getVerlet(r_cut);

	t.stop();

	return t.getwct();
}



/*! \brief Benchmark particles' forces time
 *
 * \param NN Cell list hilbert
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> double benchmark_calc_forces_hilb(T & NN, V & vd, float r_cut)
{
	//Forces timer
	timer t;
	t.start();

	calc_forces_hilb<dim,prp>(NN,vd,r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Benchmark celllist hilbert getting time
 *
 * \param NN Cell list hilbert
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 *
 * \return real time
 */
template<typename T, typename V> double benchmark_get_celllist_hilb(T & NN, V & vd, float r_cut)
{
	// Cell list timer
	timer t;
	t.start();

	//get cell list
	NN = vd.getCellList_hilb(r_cut);

	t.stop();

	return t.getwct();
}

/*! \brief Move particles in random direction but the same distance
 *
 * \param vd Distributed vector
 * \param dist Distance for which the particles are moved (note that the direction is random)
 *
 */
template<unsigned int dim, typename v_dist> void move_particles(v_dist & vd, double dist)
{
	double phi;
	double theta;
	const double pi = 3.14159;

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		phi = rand()/double(RAND_MAX);
		theta = rand()/double(RAND_MAX);
		phi *= 2.0*pi;
		theta *= pi;

		auto key = it.get();

		if(dim == 1)
			vd.getPos(key)[0] += dist*sin(theta)*cos(phi);
		else if(dim == 2)
		{
			vd.getPos(key)[0] += dist*sin(theta)*cos(phi);
			vd.getPos(key)[1] += dist*sin(theta)*sin(phi);
		}
		else if(dim == 3)
		{
			vd.getPos(key)[0] += dist*sin(theta)*cos(phi);
			vd.getPos(key)[1] += dist*sin(theta)*sin(phi);
			vd.getPos(key)[2] += dist*cos(phi);
		}

		++it;
	}
}

///////////////////////////// CONSTRUCT GRAPH //////////////////////////////

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
							  std::string y_string)
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

	if (names.size() != yp_mean.size())
		std::cerr << __FILE__ << ":" << __LINE__ << ", Error names.size() != yp_mean.size() " << names.size() << " != " << yp_mean.size() << std::endl;

	if (names.size() == 0)
		return;

	for (size_t i = 0 ; i < names.size() ; i++)
		yn2.add(names.get(i));

	for (size_t i = 0; i < xp.size() ; i++)
		x.add(xp.get(i));

	// We are assuming that yp_mean.get(g).size() are equal for each g
	y2.resize(yp_mean.get(0).size());
	y2_dev.resize(yp_mean.get(0).size());
	for (size_t r = 0; r < yp_mean.get(0).size() ; r++)
	{
		y2.get(r).resize(yp_mean.get(0).get(r).size());
		y2_dev.get(r).resize(yp_mean.get(0).get(r).size());
		for (size_t k = 0; k < yp_mean.get(0).get(r).size(); k++)
		{
			// Number of graph points
			for (size_t g = 0 ; g < yp_mean.size() ; g++)
			{
				// Put a total time
				y2.get(r).get(k).add(yp_mean.get(g).get(r).get(k));
				y2.get(r).get(k).add(yp_mean.get(g).get(r).get(k));

				y2_dev.get(r).get(k).add(yp_dev.get(g).get(r).get(k));
				y2_dev.get(r).get(k).add(yp_dev.get(g).get(r).get(k));
			}
		}
	}

	y2.save(file_mean_save);
	y2_dev.save(file_var_save);

	if (y_ref_mean.size() != 0)
	{
		// We reconstruct y and yn

		y2.clear();
		yn2.clear();

		for (size_t i = 0 ; i < yp_mean.size() ; i++)
		{
			yn2.add(names.get(i));
			yn2.add("interval");
			yn2.add("interval");
		}


		y2.resize(yp_mean.get(0).size());
		for (size_t r = 0; r < yp_mean.get(0).size(); r++)
		{
			int warning_level = -1;

			y2.get(r).resize(yp_mean.get(0).get(r).size());
			for (size_t k = 0; k < yp_mean.get(0).get(r).size(); k++)
			{

				// Number of graph points
				for (size_t g = 0 ; g < yp_mean.size() ; g++)
				{
					// Time for construction hilbert and random
					y2.get(r).get(k).add(yp_mean.get(g).get(r).get(k));
					y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(0) - 3.0*y_ref_dev.get(r).get(k).get(g));
					y2.get(r).get(k).add(y_ref_mean.get(r).get(k).get(0) + 3.0*y_ref_dev.get(r).get(k).get(g));

					warning_set(warning_level,yp_mean.get(g).get(r).get(k),y_ref_mean.get(r).get(k).get(g),y_ref_dev.get(r).get(k).get(g));
				}
			}
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

	std::string str2("<h2>2) Time to create the cell-list</h2>");
	cg.addHTML(str2);

	for (size_t i = 0; i < yp_mean.get(0).size() ; i++)
	{
		std::string chart_area;
		if (warning_vlevel.size() != 0)
			addchartarea(chart_area,warning_vlevel.get(i));
		options2.more = GC_Y_LOG + "," + GC_ZOOM + chart_area;

		options2.title = gnames.get(i);
		cg.AddLinesGraph(x,y2.get(i),yn2,options2);
	}
}

#endif /* SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_ */
