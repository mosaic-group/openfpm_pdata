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

void addUpdtateTime(GoogleChart & cg);

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
extern void StandardPerformanceGraph(std::string file_mean,
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
							  bool use_log);


#endif /* SRC_VECTOR_VECTOR_DIST_PERFORMANCE_UTIL_HPP_ */
