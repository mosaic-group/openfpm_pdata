/*
 * performance.hpp
 *
 *  Created on: Mar 9, 2016
 *      Author: yaroslav
 */

#ifndef SRC_PDATA_PERFORMANCE_CPP_
#define SRC_PDATA_PERFORMANCE_CPP_

#include <iostream>
#include <mpi.h>
#include "config.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Point_test.hpp"
#include <sstream>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

extern const char * test_dir;

// XML report
boost::property_tree::ptree pt;

#ifdef PERFORMANCE_TEST

BOOST_AUTO_TEST_SUITE( performance )

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

#include "Vector/performance/vector_dist_performance_util.hpp"
#include "Vector/performance/verlet_performance_tests.hpp"
#include "Vector/performance/cell_list_part_reorder.hpp"
#include "Vector/performance/cell_list_comp_reorder.hpp"


BOOST_AUTO_TEST_SUITE_END()

#endif

#endif /* SRC_PDATA_PERFORMANCE_CPP_ */
